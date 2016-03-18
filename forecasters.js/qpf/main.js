/* jshint node: true, esnext: true */

"use strict";

var fs = require('fs-extra');
var path = require('path');
var cp = require('child_process');
var os = require('os');

var username = require('username');
var parseArgs = require('minimist');
var debug = require('debug')('forecaster');
var pgp = require('pg-promise')({
  // Initialization Options
});

var templates = require('./config.js');
var templates = require('./templates.js');

// Parse the command line arguments
var argv = parseArgs(process.argv.slice(2), {
  string: ['output'],
  boolean: 'qsub',
  alias: {
    o: 'output'
  },
  default: {
    output: 'out'
  }
});

// Set the ouput directory
var outputDir = argv.output;

// Get the latest available state
function getLatestState() {
  var re = /^state_(\d+).rec$/;
  var dates = fs.readdirSync(outputDir)
    .filter(function (filename) {
      var stat = fs.statSync(path.join(outputDir, filename));
      return (stat.isFile() && re.test(filename));
    })
    .map(function (basename) {
      var m = re.exec(basename);
      if (m !== null) {
        return parseInt(m[1]);
      }
      return null;
    });

  return Math.max.apply(Math, dates);
}

// Read the config file
try {
  var cfg = JSON.parse(fs.readFileSync('.forecaster', 'utf8'));
} catch (ex) {
  debug('Missing or invalid config file found, generating one');
  debug(ex.message);

  // Get current time rounded to the nearest 5 minutes
  var now = Math.floor(Date.now() / 1000);
  now = now - (now % 5);

  var cfg = {
    username: username.sync(),
    obsTime: getLatestState(),
    qpfTime: now
  };
  fs.writeFileSync('.forecaster', JSON.stringify(cfg), 'utf8');
}

// Create output dir if not exists
fs.mkdirsSync(outputDir);

function renderAll(context) {
  // Render all templates
  Object.keys(templates).forEach(function (key, index) {

    var file = path.join(outputDir, key);
    debug('Rendering ' + file);

    var content = templates[key](context);
    fs.writeFileSync(file, content);
  });
}

function render(template, out, context) {
  var file = path.join(outputDir, out);
  debug('Rendering ' + file);

  var content = template(context);
  fs.writeFileSync(file, content);
}

// Generate the storm files
function generateStormfile(filePath, begin, rows) {
  //Maps the link
  var links = {};

  var length = rows.length;
  for (var i = 0; i < length; i++) {
    if (typeof links[rows[i].link_id] === 'undefined') {
      links[rows[i].link_id] = {};
      links[rows[i].link_id][begin] = 0.0;
    }

    if (rows[i].unix_time !== null) {
      links[rows[i].link_id][rows[i].unix_time] = rows[i].rain_intens;
    }
  }

  // Generate the storm files
  var outFile = fs.createWriteStream(filePath);
  outFile.write(Object.keys(links).length + os.EOL);

  Object.keys(links).forEach(function (key, index) {
    var link = links[key];
    outFile.write(key + ' ' + Object.keys(link).length + os.EOL);
    Object.keys(link).forEach(function (time, index) {
      outFile.write(time + ' ' + link[time] + os.EOL);
    });
  });

  outFile.end();
}

// Creating a new database instance from the connection details
Promise
  .all([
    pgp.connect(obsConnString),
    pgp.connect(qpfConnString)
  ])
  .then(function (clients) {
    var dbObs = clients[0];
    var dbQpf = clients[1];

    // Query the DB to get the latest OBS and QPF timestamp
    Promise
      .all([
        dbObs.one('SELECT unix_time FROM rain_maps5_index ORDER BY unix_time DESC LIMIT 1'),
        dbQpf.one('SELECT time_utc FROM hrrr_index_ldm ORDER BY time_utc DESC LIMIT 1')
      ])
      .then(function (rows) {
        var obsTime = rows[0].unix_time;
        var qpfTime = rows[1].time_utc;

        debug('lastest rainfall OBS timestamp ' + obsTime);
        debug('lastest rainfall QPF timestamp ' + qpfTime);

        // If we have new obs or qpf
        if ((cfg.obsTime < obsTime) || (cfg.qpfTime < qpfTime)) {
          var contexts = {
            obs: {
              begin: cfg.obsTime,
              end: obsTime,
              user: cfg.username,
              rainFileType: 1,
              rainFile: 'forcing_rain_obs.str',
              iniStateFile: 'state_' + cfg.obsTime + '.rec',
              endStateFile: 'state_' + obsTime + '.rec',
            },
            qpf: {
              begin: obsTime,
              end: qpfTime,
              user: cfg.username,
              rainFileType: 1,
              rainFile: 'forcing_rain_qpf.str',
              iniStateFile: 'state_' + obsTime + '.rec',
              endStateFile: 'forecast_qpf_' + qpfTime + '.rec'
            },
            s2mm: {
              begin: obsTime,
              end: obsTime + 14400 * 60,
              user: cfg.username,
              rainFileType: 4,
              rainFile: 'forcing_rain_2inches24hour.ustr',
              iniStateFile: 'state_' + obsTime + '.rec',
              endStateFile: 'forecast_2mm_' + qpfTime + '.rec'
            }
          };

          //Compute duration
          Object.keys(contexts).forEach(function (key) {
            contexts[key].duration = (contexts[key].end - contexts[key].begin) / 60.0;
          });

          // Render a new set of config files
          render(templates.gbl, 'obs.gbl', contexts.obs);
          render(templates.gbl, 'qpf.gbl', contexts.qpf);
          render(templates.gbl, 's2mm.gbl', contexts.s2mm);
          render(templates.job, 'run.job', {
            globalFiles: ['obs.gbl', 'qpf.gbl', 's2mm.gbl'],
            workingDir: path.resolve(outputDir)
          });

          debug('select OBS from ' + contexts.obs.begin + ' to ' + contexts.obs.end);
          var obsQuery = `
        SELECT unix_time, rain_intens, link_id
        FROM
          (SELECT link_id FROM materialized_env_master_km) links LEFT JOIN
          (SELECT * FROM link_rain5 WHERE unix_time >= $1 AND unix_time < $2 AND rain_intens > 0.0) rain USING (link_id)
        ORDER BY link_id, unix_time`;
          debug(obsQuery);

          debug('select QPF from ' + contexts.qpf.begin + ' to ' + contexts.qpf.end);
          var qpfQuery = `
        SELECT unix_time, rain_intens, link_id
        FROM
          (SELECT link_id FROM materialized_env_master_km) links LEFT JOIN
          (SELECT * FROM link_rain WHERE unix_time >= $1 AND unix_time < $2 AND rain_intens > 0.0) rain USING (link_id)
        ORDER BY link_id, unix_time`;
          debug(qpfQuery);

          // Generate forcing files
          Promise.all([
        dbObs.any(obsQuery, [contexts.obs.begin, contexts.obs.end])
          .then(function (rows) {
              debug('got ' + rows.length + ' OBS rainfall rows');
              generateStormfile(path.join(outputDir, contexts.obs.rainFile), contexts.obs.begin, rows);
              //dbObs.done();
            }),
        dbQpf.any(qpfQuery, [contexts.qpf.begin, contexts.qpf.end])
          .then(function (rows) {
              debug('got ' + rows.length + ' QPF rainfall rows');
              generateStormfile(path.join(outputDir, contexts.qpf.rainFile), contexts.qpf.begin, rows);
              //dbQpf.done();
            })
        ]).then(function () {

            if (!argv.qsub) return;

            // Check whether a job is already running
            if (cfg && cfg.jobId) {
              try {
                var stat = cp.execSync('qstat -j ' + cfg.jobId).toString();
                debug(stat);
              } catch (ex) {
                debug(ex);
              }
            }

            // Run the simulations
            debug('Run the simulations');
            try {
              var re = /Your job (\d*)/;

              // Queue the job
              var sub = cp.execSync('qsub ' + path.join(outputDir, 'run.job')).toString();
              var m = re.exec(sub);
              if (m !== null) cfg.jobId = parseInt(m[1]);
              debug(sub);

              // Update config file for next run
              cfg.obsTime = obsTime;
              cfg.qpfTime = qpfTime;

              fs.writeFileSync('.forecaster', JSON.stringify(cfg), 'utf8');
            } catch (ex) {
              console.error(ex.message);
            }

          }).catch(function (err) {
            return console.error('error running query', err);
          });
        }

      })
      .catch(function (err) {
        return console.error('error running query', err);
      });
  }).catch(function (err) {
    debug('ERROR', err);
  });