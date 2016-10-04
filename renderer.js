"use strict";
// This file is required by the index.html file and will
// be executed in the renderer process for that window.
// All of the Node.js APIs are available in this process.

let geneysis = require('./js/geneysis.js')

const ipc = require('electron').ipcRenderer
const dialog = require('electron').dialog
const BrowserWindow = require('electron').BrowserWindow
const shell = require('shelljs');

var app = angular.module('Geneysis', [require('angular-route'),require('angular-animate'), require('angular-ui-bootstrap'), 'electangular']);

app.config(function($routeProvider, $locationProvider) {
  $routeProvider
    .when('/', {
        template: '<div></div>',
        controller: 'StartController',
    })
    .when('/import', {
    templateUrl: './templates/import.html',
    controller: 'ImportController',
  })
});
app.service("globalSettings", [function() {
  var observerCallbacks = [];

  //register an observer
  var registerObserverCallback = function(callback){
    observerCallbacks.push(callback);
  };

  //call this when you know 'foo' has been changed
  var notifyObservers = function(){
    angular.forEach(observerCallbacks, function(callback){
      callback();
    });
  };
  
  
  return {
    project: undefined,
    config: undefined,
    observerCallbacks: observerCallbacks,
    registerObserverCallback: registerObserverCallback,
    notifyObservers: notifyObservers
  }
}]);

app.controller("MainController", ["$scope","globalSettings","$window", function($scope, globalSettings, $window) {
  $scope.sidebarOpen = true;
  $scope.mainWidth = 0;
  $scope.project = globalSettings.project;
  $scope.setMainWidth = function() {
    if ($scope.sidebarOpen) {
      $scope.mainWidth = $window.innerWidth - 304;
    }
    else {
      $scope.mainWidth = $window.innerWidth;
    }
  }
  $scope.setMainWidth();
  angular.element($window).bind('resize',function() {
    $scope.setMainWidth();
    $scope.$digest();
  })
  
  
  globalSettings.registerObserverCallback(function() {
    console.log(globalSettings);
    $scope.project = globalSettings.project;
  })
  
}]);

app.controller("StartController", ["$rootScope","$scope","$window","ipc","globalSettings", function($rootScope,$scope, $window,ipc,globalSettings) {
  
  ipc.send({task:"get_current_project"}); // ask for the current project
  $rootScope.$on('electron-msg', (event, msg) => { // listen for response
    if (msg.task == "load_project") { // make sure the message is really for you
      console.log(msg)
      globalSettings.project = msg.project;
      globalSettings.project.created = new Date(globalSettings.project.created);
      globalSettings.project.updated = new Date(globalSettings.project.updated);
      if (globalSettings.project.events.length == 0) {
        $window.location.hash = '#import';
//        console.log($window.location);
      }
      globalSettings.notifyObservers();
    }
    
  });
  
}]);

app.controller("ImportController", ["$rootScope","$scope","electron","ipc","globalSettings","$location", function($rootScope, $scope, electron, ipc, globalSettings,$location) {
  if (globalSettings.project == undefined) {
    $location.path("/");
  }
  else {
    $scope.project = globalSettings.project;
  }
  console.log("STARTING!");
  $rootScope.$on('electron-msg', (event, msg) => {
    
    switch(msg.task) {
        case "load_genbank":
          var files = msg.files;
          for (let i=0;i<files.length;i++) {
            geneysis.loadPhage($scope.project.path, files[i], function(results, err) {
              if (err) {
                console.log(err)
              }
              try {
                globalSettings.project = JSON.parse(results);
                $scope.project = globalSettings.project;
                globalSettings.notifyObservers();
                console.log($scope.project)
              }
              catch(error) {
                console.error(error);
              }
            });
          }
          break;
        
        default:
          console.log(msg);
          break;
    }
    
  });
    
  $scope.openFiles = function() {
   ipc.send({task:"get_genbank"});
  }
  
  
}]);
