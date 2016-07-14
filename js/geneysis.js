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

var Geneysis = {
  
  createNewDir : function(wd, callbackFn) {
    var options = {
      mode: 'text',
      scriptPath: './python/',
      args: ['--wd', wd, 'create']
    };
    PythonShell.run(scriptPath, options, function(err, results) {
      callbackFn(results, err)
    })  
  },
  
  load : function(wd, file, callbackFn) {
    var options = {
      mode: 'text',
      scriptPath: './python/',
      args: ['--wd', wd, 'import', '--file', file]
    };
    PythonShell.run(scriptPath, options, function(err, results) {
      callbackFn(results, err)
    })  
  },
  
 
}

 module.exports = Geneysis;