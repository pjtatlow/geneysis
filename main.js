const electron = require('electron');
// Module to control application life.
const app = electron.app;
// Module to create native browser window.
const BrowserWindow = electron.BrowserWindow;
// Module to communicate with renderer process.
const ipcMain = require('electron').ipcMain;
// Module to open diaglogs
const dialog = require('electron').dialog;
// Modules to allow access to shell commands
const shell = require('shelljs');
const exec = require('child_process').exec;

// Keep a global reference of the window object, if you don't, the window will
// be closed automatically when the JavaScript object is garbage collected.
let mainWindow;
let projectsWindow;
let project = undefined;

function openSelectProject () {
  
  // Create the browser window.
  projectsWindow = new BrowserWindow({width: 400, height: 600, resizable:false});

  projectsWindow.loadURL(`file://${__dirname}/selectProject.html`);

  // Emitted when the window is closed.
  projectsWindow.on('closed', function () {
    // Dereference the window object, usually you would store windows
    // in an array if your app supports multi windows, this is the time
    // when you should delete the corresponding element.
    projectsWindow = null;
  });
  
}

function openGeneysis() {

  // Create the browser window.
  mainWindow = new BrowserWindow({width: 1400, height: 900});

  mainWindow.loadURL(`file://${__dirname}/index.html`);

  // Open the DevTools.
  mainWindow.webContents.openDevTools();

  // Emitted when the window is closed.
  mainWindow.on('closed', function () {
    // Dereference the window object, usually you would store windows
    // in an array if your app supports multi windows, this is the time
    // when you should delete the corresponding element.
    mainWindow = null;
  });
  
  if (projectsWindow != null) {
    projectsWindow.close();
  }  
  
}

function startGeneysis () {

  if (project == undefined) {
    openSelectProject();
  }
  else {
    openGeneysis();
  }
  //check for all dependencies
  var missing = "";
  if (!shell.which("clustalo")) {
    missing += "Clustal Omega (clustalo): http://clustal.org/omega/\n";
  }
  if (!(shell.which("blastp") && shell.which("makeblastdb"))) {
    missing += "Blast: http://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local/\n"
  }
  if (!shell.which("python")) {
    missing += "Python 2.7+: https://www.python.org/\n";
  }
  exec('python --version', (error, stdout, stderr) => { // check python version 2.7+
    if (error) {
      console.error("FAILED TO CHECK PYTHON VERSION");
      app.quit();
    }
    stderr = stderr.match(/Python (\d\.\d).+/)[1]
    if (stderr != '2.7') {
      missing += "Python 2.7+: https://www.python.org/\n";
    }

    exec('python -c "import sqlite3"; echo $?', (error, stdout, stderr) => { // make sure sqlite3 module is installed
      if (error) {
          dialog.showMessageBox(mainWindow, {type:"error",title:"ERROR",message:'python -c "import sqlite3"; echo $?',buttons:["OK"]}, function() {
            mainWindow.close();
          });
      }        
      else {
        stdout = stdout.replace("\n","")
        if (stdout != "0") {
          missing += "Missing SQLite3 Module: try \"pip install pysqlite\"\n";
        }
      }

      exec('python -c "from Bio import SeqIO"; echo $?', (error, stdout, stderr) => { // make sure sqlite3 module is installed
        if (error) {
          dialog.showMessageBox(mainWindow, {type:"error",title:"ERROR",message:'Failed to run python -c "from Bio import SeqO"; echo $?',buttons:["OK"]}, function() {
            mainWindow.close();
          });
        }        
        else {
          stdout = stdout.replace("\n","")
          if (stdout != "0") {
            missing += "Missing BioPython Module: try \"pip install biopython\"\n";
          }
        }
        if (missing.length > 0) {
          dialog.showMessageBox(mainWindow, {type:"error",title:"Missing Dependencies",message:"Please install the following dependencies:",detail:missing,buttons:["OK"]}, function() {
            mainWindow.close();
          });
        }
      });
    });




  });


}

// This method will be called when Electron has finished
// initialization and is ready to create browser windows.
// Some APIs can only be used after this event occurs.
app.on('ready', startGeneysis);

// Quit when all windows are closed.
app.on('window-all-closed', function () {
  app.quit();
});

// In this file you can include the rest of your app's specific main process
// code. You can also put them in separate files and require them here.

ipcMain.on('electron-msg', (event, msg) => {
  //handle incoming message here
  console.log(msg);
  if (msg.task == "get_genbank") { // open dialog to allow user to select files
    
    let options = { 
      properties: ['openFile', 'multiSelections'],
      filters: [
        {name: 'Genbank', extensions: ['.gbk', '.gb', '.genbank']},
        {name: 'Others', extensions: ['*']}
      ],
    };
    
    if (process.platform !== 'darwin') {
      dialog.showOpenDialog(options, function (selected_files) {
        if (selected_files)  {
          event.sender.send('electron-msg', {
            task: "load_genbank",
            files: selected_files
          })
        }
      });

    }
    else { // use sheet dialog
      const window = BrowserWindow.fromWebContents(event.sender)
      const files = dialog.showOpenDialog(window, options, function (selected_files) {
        if (selected_files)  {
          event.sender.send('electron-msg', {
            task: "load_genbank",
            files: selected_files
          })
        }
      });
    }

    
  }
  else if (msg.task == "get_projects") {
    shell.cd();
    var home = shell.pwd().stdout;
    shell.cd('.geneysis');
    if (shell.pwd().stdout == home) {
      shell.mkdir('-p',[".geneysis/projects",".geneysis/config"]);
      shell.cd('.geneysis');
    }
    var wd = shell.pwd().stdout
    var results = shell.ls('-d','projects/*');
    var projects = []
    for (let i = 0; i < results.length; i++) {
      let obj = {};
      obj.name = results[i].substring(results[i].lastIndexOf('/') + 1);
      obj.path = wd + "/" + results[i] + "/";
      projects.push(obj);
    }
    
    
    event.sender.send('electron-msg', {
      task: "projects",
      projects: projects,
    })    
  }
  else if (msg.task == "open_project") {
    project = msg.project;
    openGeneysis();
  }
  else if (msg.task == "get_current_project") {
    event.sender.send('electron-msg', {
      task: "load_project",
      project: project
    });
  }
});  