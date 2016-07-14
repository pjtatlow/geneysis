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
        controller: 'MainController',
    })
    .when('/home', {
    templateUrl: './templates/home.html',
    controller: 'HomeController',
  })
});
app.service("settings", [function() {
  return {
    wd: undefined,
    project: undefined,
    
  }
}]);

app.controller("MainController", ["$rootScope","$scope","$window","ipc","settings", function($rootScope,$scope, $window,ipc,settings) {
  ipc.send({task:"get_current_project"});
  $rootScope.$on('electron-msg', (event, msg) => {
    if (msg.task == "load_project") {
      console.log(msg);
    }
    
  });
  
}]);

app.controller("HomeController", ["$rootScope","$scope","electron","ipc", function($rootScope, $scope, electron, ipc) {
  
  $rootScope.$on('electron-msg', (event, msg) => {
    
    switch(msg.task) {
        case "load_genbank":
          
          break;
        
        default:
          console.log(msg);
          break;
    }
    
  });
  
  $scope.name = "YAY!"
  
  $scope.files = []
  
  $scope.$watch("files", function() {
    console.log($scope.files)
  })
  
  $scope.newDir = function() {
    geneysis.createNewDir("/home/pjtatlow/.geneysis/");
  }
  
  $scope.openFiles = function() {
   ipc.send({task:"get_genbank"});
  }
  
  
}]);
