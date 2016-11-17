let PythonShell = require('python-shell');
const scriptPath = "geneysis.py";
var requireKeys = function(obj, keys) {
  for (let i = 0; i < keys.length; i++) {
    if (obj[keys[i]] === undefined || obj[keys[i]] === null) {
      throw {"error":"MISSING KEYS", "given":obj, "required_keys":keys}
    }
  }
  return true;
};

// need more functions here
var Geneysis = {
  
  createNewProject : function(wd, callbackFn) {
    console.log(wd);
    var options = {
      mode: 'text',
      scriptPath: `${__dirname}/../python/`,
      args: ['create', '--wd', wd]
    };
    PythonShell.run(scriptPath, options, function(err, results) {
      callbackFn(results, err)
    })  
  },
  
  loadPhage : function(wd, file, callbackFn) {
    var options = {
      mode: 'text',
      scriptPath: `${__dirname}/../python/`,
      args: ['import', '--wd', wd, '--file', file]
    };
    PythonShell.run(scriptPath, options, function(err, results) {
      callbackFn(results, err)
    })  
  },
  
 
}

 module.exports = Geneysis;
