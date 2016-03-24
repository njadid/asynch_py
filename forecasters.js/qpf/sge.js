/* jshint node: true, esnext: true */

"use strict";

var cp = require('child_process');
var debug = require('debug')('sge');
var parseString = require('xml2js').parseString;

module.exports = {

  // Submit a new job
  qsub: function(jobFile, depJobNames) {

    return new Promise(function (fulfill, reject) {
      var re = /Your job (\d*)/;

      //Format the dependency options
      var cmd = 'qsub ';
      if (depJobNames) {
        cmd += '-hold_jid ' + depJobNames.join(',') + ' ';
      }
      cmd += jobFile;

      debug(cmd);
      cp.exec(cmd, function (err, stdout, stderr) {
        if (err)  {
          reject(err);
        } else {
          debug(stdout);
          debug(stderr);
          var m = re.exec(stdout);
          if (m !== null) {
            fulfill(parseInt(m[1]));
          } else {
            reject(new Error('Failed to parse the qsub response'));
          }
        }
      });
    });
  },

  // Check job status
  qstat: function(jobName) {
    return new Promise(function (fulfill, reject) {
      var cmd = 'qstat -xml -u $USER';

      debug(cmd);
      var stat = cp.exec(cmd, /*{encoding: 'utf16'},*/ function (err, stdout, stderr) {
        if (err)  {
          reject(err);
        } else {
          debug(stdout);
          debug(stderr);
          parseString(stdout, {trim: true, emptyTag: null, explicitArray: false}, function (err, result) {
            if (err)  {
              reject(err);
            }
            if (result.job_info && result.job_info.queue_info && result.job_info.queue_info.job_list) {
              fulfill(result.job_info.queue_info.job_list.find(function (job) { return job.JB_name === jobName; }));
            } else {
              fulfill(undefined);
            }
          });
        }
      });
    });
  },
  
  // Delete a job
  qdel: function(jobName) {
    return new Promise(function (fulfill, reject) {
      var cmd = 'qdel ' + jobName;

      debug(cmd);
      var stat = cp.exec(cmd, /*{encoding: 'utf16'},*/ function (err, stdout, stderr) {
        if (err)  {
          reject(err);
        } else {
          fulfill();
        }
      });
    });
  }
};
