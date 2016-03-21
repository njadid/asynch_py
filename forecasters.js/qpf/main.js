/* jshint node: true, esnext: true */

"use strict";

var fs = require('fs-extra');
var path = require('path');
var cp = require('child_process');
var os = require('os');

var username = require('username');
var parseArgs = require('minimist');
var debug = require('debug')('forecaster');
var pgp = require('pg-promise-strict');

var config = require('./config.js');
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

// For each row add intensity to the map of rain data
function mapRainLink(begin, links) {

  return function(row) {
    if (typeof links[row.link_id] === 'undefined') {
      links[row.link_id] = {};
      links[row.link_id][begin] = 0.0;
    }

    if (row.unix_time !== null) {
      links[row.link_id][row.unix_time] = row.rain_intens;
    }
  };
}

// Generate the storm files
function generateStormfile(filePath, begin, links) {
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
    pgp.connect(config.obsConn),
    pgp.connect(config.qpfConn)
  ])
  .then(function (clients) {
    var obsClient = clients[0];
    var qpfClient = clients[1];

    // Query the DB to get the latest OBS and QPF timestamp
    Promise
      .all([
        obsClient.query('SELECT unix_time FROM rain_maps5_index ORDER BY unix_time DESC LIMIT 1').fetchUniqueValue(),
        qpfClient.query('SELECT time_utc FROM hrrr_index_ldm ORDER BY time_utc DESC LIMIT 1').fetchUniqueValue()
      ])
      .then(function (results) {
        var obsTime = results[0].value;
        var qpfTime = results[1].value;

        debug('lastest rainfall OBS timestamp ' + obsTime);
        debug('lastest rainfall QPF timestamp ' + qpfTime);

        // If we have new obs of qpf
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
              outHydrographsFile: {
                file: 'hydrographs.dat'
              },
              outPeaksFile: {
                file: 'peaks.pea'
              }
            },
            qpf: {
              begin: obsTime,
              end: qpfTime,
              user: cfg.username,
              rainFileType: 1,
              rainFile: 'forcing_rain_qpf.str',
              iniStateFile: 'state_' + obsTime + '.rec',
              endStateFile: 'forecast_qpf_' + qpfTime + '.rec',
              outHydrographsDb: {
                file: 'save_hydrograph.dbc',
                table: 'qpf.hydrographs'
              },
              outPeaksDb: {
                file: 'save_peaks.dbc',
                table: 'qpf.peaks'
              }
            },
            s2in: {
              begin: obsTime,
              end: obsTime + 14400 * 60,
              user: cfg.username,
              rainFileType: 4,
              rainFile: 'forcing_rain_2inches24hour.ustr',
              iniStateFile: 'state_' + obsTime + '.rec',
              endStateFile: 'forecast_2in_' + qpfTime + '.rec',
              outHydrographsDb: {
                file: 'save_hydrograph.dbc',
                table: 'whatif_2in24h.hydrographs'
              },
              outPeaksDb: {
                file: 'save_peaks.dbc',
                table: 'whatif_2in24h.peaks'                
              }
            }
          };

          //Compute duration
          Object.keys(contexts).forEach(function (key) {
            contexts[key].duration = (contexts[key].end - contexts[key].begin) / 60.0;
          });

          // Render a new set of config files
          render(templates.gbl, 'obs.gbl', contexts.obs);
          render(templates.gbl, 'qpf.gbl', contexts.qpf);
          render(templates.gbl, 's2in.gbl', contexts.s2in);
          render(templates.job, 'run.job', {
            globalFiles: ['obs.gbl', 'qpf.gbl', 's2in.gbl'],
            workingDir: path.resolve(outputDir)
          });

          debug('select OBS from ' + contexts.obs.begin + ' to ' + contexts.obs.end);
          var obsQuery = `
SELECT unix_time, rain_intens, link_id
FROM
  (SELECT link_id FROM materialized_env_master_km) links LEFT JOIN
  (SELECT * FROM link_rain5 WHERE unix_time >= $1 AND unix_time < $2 AND rain_intens > 0.0) rain USING (link_id)
ORDER BY link_id, unix_time`;
          //debug(obsQuery);

          debug('select QPF from ' + contexts.qpf.begin + ' to ' + contexts.qpf.end);
          var qpfQuery = `
SELECT unix_time, rain_intens, link_id
FROM
  (SELECT link_id FROM materialized_env_master_km) links LEFT JOIN
  (SELECT * FROM link_rain WHERE unix_time >= $1 AND unix_time < $2 AND rain_intens > 0.0) rain USING (link_id)
ORDER BY link_id, unix_time`;
          //debug(qpfQuery);

          var obsLinks = {}, qpfLinks = {};
          Promise
            .all([
              obsClient.query(obsQuery, [contexts.obs.begin, contexts.obs.end])
                .onRow(mapRainLink(contexts.obs.begin, obsLinks)),
              qpfClient.query(qpfQuery, [contexts.qpf.begin, contexts.qpf.end])
                .onRow(mapRainLink(contexts.qpf.begin, qpfLinks))
            ])
            .then(function (results) {

              debug('got ' + results[0].rowCount + ' OBS rainfall rows');
              debug('got ' + results[1].rowCount + ' QPF rainfall rows');

              generateStormfile(path.join(outputDir, 'forcing_rain_obs.str'), contexts.obs.begin, obsLinks);
              generateStormfile(path.join(outputDir, 'forcing_rain_qpf.str'), contexts.qpf.begin, qpfLinks);

              obsClient.done();
              qpfClient.done();

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

            })
            .catch(function (err) {
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
