/* jshint node: true, esnext: true */

"use strict";

var fs = require('fs-extra');
var path = require('path');
var os = require('os');

var username = require('username');
var parseArgs = require('minimist');
var debug = require('debug')('forecaster');
var pgp = require('pg-promise-strict');

var config = require('./config.js');
var templates = require('./templates.js');
var sge = require('./sge.js');
var stormfile = require('./stormfile.js');
//var scenarii = require('./scenarii.js');

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

// Create output dir if not exists
fs.mkdirsSync(outputDir);

//Render a template
function render(template, out, context) {
  var file = path.join(outputDir, out);
  debug('Rendering ' + file);

  var content = template(context);
  fs.writeFileSync(file, content);
}

//Check whether a simulation is already running
sge.qstat('IFC_QPE')
.then(function (stat) {
  debug(stat);
})
.then(function (stat) {
  // Creating a new database instance from the connection details
  return pgp.connect(config.obsConn).then(function (client) {
    // Query the DB to get the latest obs timestamp
    return client.query('SELECT unix_time FROM rain_maps5_index ORDER BY unix_time DESC LIMIT 1');
  });
})
.then(function (query){
  return query.fetchUniqueValue();
}).then(function (result) {
  var currentTime = getLatestState(),
    obsTime = result.value,
    user = username.sync();

  debug('current timestamp ' + currentTime);
  debug('lastest rainfall QPE timestamp ' + obsTime);

  // If we have new obs ready
  if (currentTime >= obsTime)
    throw 'Simulation is already up to date';

  var context = {
    user: user,
    begin: currentTime,
    end: obsTime,
    duration: (obsTime - currentTime) / 60,
    iniStateFile: 'state_' + currentTime + '.rec',
    endStateFile: 'state_' + obsTime + '.rec',
    rainFileType: 1,
    rainFile: 'forcing_rain_obs.str',
    outHydrographsFile: {
      file: 'hydrographs.dat'
    },
    outPeaksFile: {
      file: 'peaks.pea'
    }
  };

  // Render a new set of config files
  render(templates.gbl, 'qpe.gbl', context);
  render(templates.job, 'qpe.job', {
    name: 'IFC_QPE',
    globalFile: 'qpe.gbl',
    workingDir: path.resolve(outputDir)
  });

  // Get the QPE data and generate the stormfile
  var links = {};
  return result.client.query(config.obsQuery, [context.begin])
    .onRow(stormfile.mapLink(context.begin, links))
  .then(function (result) {
    debug('got ' + result.rowCount + ' QPE rainfall rows');
    result.client.done();

    stormfile.generate(path.join(outputDir, context.rainFile), context.begin, links);

    // Run the simulations
    debug('Queue the simulations');
    return sge.qsub(path.join(outputDir, 'qpe.job'));
  })
  .then(function (jobId) {
    debug('Generate the QPF forecast simulations');

    links = {};
    // Creating a new database instance from the connection details
    return pgp.connect(config.qpfConn).then(function (client) {
      // Query the DB to get the latest obs timestamp
      return client.query(config.qpfQuery, [context.begin])
        .onRow(stormfile.mapLink(context.begin, links));
    }).then(function (result) {
      debug('got ' + result.rowCount + ' QPF rainfall rows');
      result.client.done();

      var context = {
        user: user,
        begin: currentTime,
        end: obsTime + 14400 * 60,
        duration: 14400,
        iniStateFile: 'state_' + obsTime + '.rec',
        endStateFile: 'forecast_qpf.rec',
        rainFileType: 1,
        rainFile: 'forcing_rain_qpf.str',
        outHydrographsDb: {
          file: 'save_hydrograph.dbc',
          table: 'qpf.hydrographs'
        },
        outPeaksDb: {
          file: 'save_peaks.dbc',
          table: 'qpf.peaks'
        }
      };

      // Render a new set of config files
      render(templates.gbl, 'qpf.gbl', context);
      render(templates.job, 'qpf.job', {
        name: 'IFC_QPF',
        globalFile: 'qpf.gbl',
        workingDir: path.resolve(outputDir)
      });

      stormfile.generate(path.join(outputDir, context.rainFile), context.begin, links);

      // Run the simulations
      debug('Queue the simulations');
      return sge.qsub(path.join(outputDir, 'qpf.job'), ['IFC_QPE']);
    });
  })
  .then(function (jobId) {
    debug('Generate the scenarii forecast simulations');

    var context = {
      user: user,
      begin: obsTime,
      end: obsTime + 14400 * 60,
      duration: 14400,
      rainFileType: 4,
      rainFile: 'forcing_rain_2inches24hour.ustr',
      iniStateFile: 'state_' + obsTime + '.rec',
      endStateFile: 'forecast_2in.rec',
      outHydrographsDb: {
        file: 'save_hydrograph.dbc',
        table: 'whatif_2in24h.hydrographs'
      },
      outPeaksDb: {
        file: 'save_peaks.dbc',
        table: 'whatif_2in24h.peaks'
      }
    };
    // Render a new set of config files
    render(templates.gbl, 'whatif_2in24h.gbl', context);
    render(templates.job, 'whatif_2in24h.job', {
      name: 'IFC_WHATIF_2IN24H',
      globalFile: 'whatif_2in24h.gbl',
      workingDir: path.resolve(outputDir)
    });

    // Run the simulations
    debug('Queue the simulations');
    return sge.qsub(path.join(outputDir, 'whatif_2in24h.job'), ['IFC_QPE']);
  });
})
.catch(function (err) {
  return console.error(err);
}).then(function(){
  process.exit();
});