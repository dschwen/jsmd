/*************************************************************
 *
 *  MathJax/extensions/TeX/bbox.js
 *  
 *  This file implements the \bbox macro, which creates an box that
 *  can be styled (for background colors, and so on).  You can include
 *  an optional dimension that tells how much extra padding to include
 *  around the bounding box for the mathematics, or a color specification 
 *  for the background color to use, or both.  E.g.,
 *  
 *    \bbox[2pt]{x+y}        %  an invisible box around x+y with 2pt of extra space
 *    \bbox[green]{x+y}      %  a green box around x+y
 *    \bbox[green,2pt]{x+y}  %  a green box with 2pt of extra space
 *  
 *  This file will be loaded automatically when \bbox is first used.
 *
 *  ---------------------------------------------------------------------
 *  
 *  Copyright (c) 2011 Design Science, Inc.
 * 
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

MathJax.Hub.Register.StartupHook("TeX Jax Ready",function () {
  var VERSION = "0.9";

  var TEX = MathJax.InputJax.TeX, MML = MathJax.ElementJax.mml;
  TEX.Definitions.macros.bbox = "BBox";
  TEX.Parse.Augment({
    BBox: function (name) {
      var bbox = this.GetBrackets(name),
          math = this.ParseArg(name);
      var parts = bbox.split(/,/), def = {};
      for (var i in parts) {
        var part = parts[i].replace(/^\s+/,'').replace(/\s+$/,'');
        var match = part.match(/^(\.\d+|\d+(\.\d*)?)(pt|em|ex|mu|px)$/);
        if (match) {
          def.height = def.depth = "+"+match[1]+match[3];
          def.width = "+"+(2*match[1])+match[3];
          def.lspace = match[1]+match[3];
        } else if (part.match(/^([a-z0-9]+|\#[0-9a-f]{6}|\#[0-9a-f]{3})$/i)) {
          def.mathbackground = part;
        } else if (part.match(/^[-a-z]+:/i)) {
          def.style = part;
        } else {
          TEX.Error("'"+part+"' doesn't look like a color or a padding dimension");
        }
      }
      this.Push(MML.mpadded(math).With(def));
    }
  });

  MathJax.Hub.Startup.signal.Post("TeX bbox Ready");

});

MathJax.Ajax.loadComplete("[MathJax]/extensions/TeX/bbox.js");
