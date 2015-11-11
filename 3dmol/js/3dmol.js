/* 3Dmol functions 
*
*
*/

   //protein object
    var protein = {
        surface: $3Dmol.SurfaceType.SAS,
        opacity: 1,
        min_isoval: -5,
        max_isoval: 5,
        colorScheme: "RWB",
        volumedata: null
    };

    var volumedata = null;
    var glviewer = null;
    var labels = [];
   

    var addLabels = function() {
        var atoms = glviewer.getModel().selectedAtoms({
            atom : "CA"
        });
        for ( var a in atoms) {
            var atom = atoms[a];

            var l = glviewer.addLabel(atom.resn + " " + atom.resi, {
                inFront : true,
                fontSize : 12,
                position : {
                    x : atom.x,
                    y : atom.y,
                    z : atom.z
                }
            });
            atom.label = l;
            labels.push(atom);
        }
    };
    
    var removetheLabels = function() {
        for (var i = 0; i < labels.length; i++) {
        var atom = labels[i]
        glviewer.removeLabel(atom.label)
        delete atom.label
        }
        //console.log(labels)
        
        labels = []

        };

    var atomcallback = function(atom, viewer) {
        if (atom.clickLabel === undefined
                || !atom.clickLabel instanceof $3Dmol.Label) {
            atom.clickLabel = viewer.addLabel(atom.elem + atom.serial, {
                fontSize : 14,
                position : {
                    x : atom.x,
                    y : atom.y,
                    z : atom.z
                },
                backgroundColor: "gray"
            });
            atom.clicked = true;
        }

        //toggle label style
        else {

            //if (atom.clicked) {
            //  var newstyle = atom.clickLabel.getStyle();
            //  newstyle.backgroundColor = 0x66ccff;

            //  viewer.setLabelStyle(atom.clickLabel, newstyle);
            //  atom.clicked = !atom.clicked;
            //}
            if (atom.clicked) {
                viewer.removeLabel(atom.clickLabel);
                delete atom.clickLabel;
                atom.clicked = false;
            }

        }
    };
    
    var glviewer;
    $(document).ready(function() {
        glviewer = $3Dmol.createViewer("gldiv", {
        defaultcolors : $3Dmol.rasmolElementColors
        });
        glviewer.setBackgroundColor("black");

    });
    
    var fileselected = function(files, func){
        
        
        readText(files, func);

        
    };
     
    var addpqr = function(data){
        
        //moldata = data = $("#moldata_pdb_large").val();
        //console.log(data); //see contents of file
        receptorModel = m = glviewer.addModel(data, "pqr");

        atoms = m.selectedAtoms({});

        for ( var i in atoms) {
            var atom = atoms[i];
            atom.clickable = true;
            atom.callback = atomcallback;
        }

        glviewer.mapAtomProperties($3Dmol.applyPartialCharges);
        glviewer.zoomTo();
        glviewer.render();
        
        };
        
        
    var addcube = function (volumedata){
        //protein.volumedata = volumedata;
        window.volumedata = new $3Dmol.VolumeData(volumedata, "cube");
        //volumedata = $("#volumetric_data").val();
        //glviewer.addIsosurface(volumedata, {isoval: -5, color:"red", smoothness: 10})
        //glviewer.addIsosurface(volumedata, {isoval: 5, color:"blue", smoothness: 1})
        
        
        glviewer.render();
        create_surface();
        };
    
    var backbone = function (){
        var atoms = glviewer.getModel().selectedAtoms({
            });
        for ( var i = 0; i < atoms.length; i++) {
            var atom = atoms[i];
        if (atom.atom == "H")
        //    delete atom
        //if (atom == "O")
        //    delete atom
        //if (atom.atom == "CA")
        atoms.splice(i,1);
        }
    }
   
    var readText = function(input,func) {
        
        if(input.length > 0) {
            var file = input[0];
            var reader = new FileReader();
            reader.onload = function(evt) {
                func(evt.target.result,file.name);
            };
            reader.readAsText(file); //needs to be type Blob
            $(input).val('');
            
        }

    };
    
    var distance = function(atom1, atom2) {
        m = glviewer.getModel(0);
        myatoms = m.selectedAtoms({});
        //console.log(myatoms)
        for ( var i in myatoms) {
        var myatom = myatoms[i];
        myatom.clickable = true;
    }   
        myatom.onclick = console.log(myatom)
    };

    /*update surface based on selected action 
    * 0 -  
    * 1 - change surface
    * 2 - set translucent
    * 3 - set opaque
    */
    function update_surface(action){
        var e = document.getElementById("selected_surface");
        var x = e.options[e.selectedIndex].value;
        glviewer.removeSurface(surf);
        switch (action){
            case 1:
                if (x == 'SAS')
                   protein.surface = $3Dmol.SurfaceType.SAS;
                else if (x == 'SES')
                    protein.surface = $3Dmol.SurfaceType.SES;
                else if (x == 'VDW')
                    protein.surface = $3Dmol.SurfaceType.VDW;
                break;
            case 2:
                protein.opacity = 0.70;
                break;
            case 3: 
                protein.opacity = 1;
                break;
            case 4:
                protein.min_isoval = -5;
                protein.max_isoval = 5;
                break;
            
            default:
                break;
        }
        
            set_color();
        }

        function set_vis(){
            var f = document.getElementById("selected_vis");
            var y = f.options[f.selectedIndex].value;
            vis=y;

            if(y=="stick"){ glviewer.setStyle({},{stick:{}}); glviewer.render();}
            if(y=="line"){glviewer.setStyle({},{line:{}}); glviewer.render();}
            if(y=="cross"){glviewer.setStyle({},{cross:{linewidth:5}}); glviewer.render();}
            if(y=="sphere"){glviewer.setStyle({},{sphere:{}}); glviewer.render();}
            if(y=="cartoon"){glviewer.setStyle({hetflag:false},{cartoon:{color: 'spectrum'}}); glviewer.render();}
        }

        function set_color(){
        //inefficient -- need to fix!
        //want to set as protein attribute
        var f = document.getElementById("selected_scheme");
        var y = f.options[f.selectedIndex].value;
        protein.colorScheme=y;
        
        if(protein.colorScheme=="RWB")
            volscheme_to_use = new $3Dmol.Gradient.RWB(protein.min_isoval,protein.max_isoval);
        else if(protein.colorScheme=="ROYGB")
            volscheme_to_use = new $3Dmol.Gradient.ROYGB(protein.min_isoval,protein.max_isoval);
        else if(protein.colorScheme=="BWR")
            volscheme_to_use = new $3Dmol.Gradient.Sinebow(protein.min_isoval,protein.max_isoval);
        
        surf = glviewer.addSurface(protein.surface, {opacity:protein.opacity, voldata: volumedata, volscheme: volscheme_to_use});
        }
        
        //starts program with SAS surface
        function create_surface(){
            volscheme_to_use = new $3Dmol.Gradient.RWB(protein.min_isoval,protein.max_isoval);
            surf = glviewer.addSurface(protein.surface, {opacity:protein.opacity, voldata: volumedata, volscheme: volscheme_to_use});
        }

        //Turn on the surface for the current selected surface
        function on_surface(){
            var e = document.getElementById("selected_surface");
            var x = e.options[e.selectedIndex].value;
            if (x == 'SAS')
                protein.surface = $3Dmol.SurfaceType.SAS;
            else if (x == 'SES')
                protein.surface = $3Dmol.SurfaceType.SES;
            else if (x == 'VDW')
                protein.surface = $3Dmol.SurfaceType.VDW;

            set_color();
        }

        //change output for min_isoval range, not perfect
        function set_min_isoval(min_val) {
            document.querySelector('#min_isoval').value = min_val;
            protein.min_isoval = min_val;
            update_surface(0);
        }

        //change output for max_isoval range, not perfect
        function set_max_isoval(max_val) {
            document.querySelector('#max_isoval').value = max_val; 
            protein.max_isoval = max_val;
            update_surface(0);
        }

        //reset min and max isovals
        //does not move slider---probably should fix this
        function reset_vals() {
            set_min_isoval2(-5);
            set_max_isoval2(5);
            document.getElementById("min_isoval2").value = "-5";
            document.getElementById("max_isoval2").value = "5";
            update_surface(0);
            return false;
        }
        
 //change output for min_isoval range, not perfect
        function set_min_isoval2(min_val) {
            document.getElementById("min_isoval").innerHTML = min_val;
            protein.min_isoval = Number(min_val);
            console.log(document.getElementById('min_isoval').value);
            update_surface(0);
        }

        //change output for max_isoval range, not perfect
        function set_max_isoval2(max_val) {
            document.getElementById("max_isoval").innerHTML = max_val;
            protein.max_isoval = Number(max_val);
            update_surface(0);
        }

        function getpqr(jobid){
            var xhr = new XMLHttpRequest();
            //jobid = 14357857643;
            //url = "@website@tmp/"+jobid+"/"+jobid+".pqr";
            url = "../3dmol/files/1fas.pqr";
            xhr.open("GET", url);
            //xhr.responseType = 'blob';

            xhr.onload = function(e) {
              if (this.status == 200) {
                // Note: .response instead of .responseText
                //var blob = new Blob([this.response], {type: 'text/plain'});
               //readText(this.response);
               addpqr(this.response);
              }
              
            };
            xhr.send(null);
            
        }

        function getcube(jobid){
            var xhr = new XMLHttpRequest();
            //xhr.open("GET", "@website@tmp/"+jobid+"/"+jobid+".cube");
            xhr.open("GET", "../3dmol/files/1fas.cube");
            //xhr.responseType = 'blob';

            xhr.onload = function(e) {
              if (this.status == 200) {
                // Note: .response instead of .responseText
                //var blob = new Blob([this.response], {type: 'text/plain'});
               //readText(this.response);
               addcube(this.response);
              }
              
            };
            xhr.send(null);
            
        }

var surfaceOn = true
function toggleSurface(){
    if(surfaceOn){
        surfaceOn = false
        on_surface()
    }
    else{
        surfaceOn = true
        glviewer.removeSurface(surf)
    }
}

var surfaceOpacity = true
function toggleOpacity(){
    if(surfaceOpacity){
        update_surface(3)
        surfaceOpacity = false
    }
    else{
        update_surface(2)
        surfaceOpacity = true
    }
}

var modelLabels = false
function toggleLabels(){
    if(modelLabels){
        removetheLabels(glviewer);
        glviewer.render();
        modelLabels = false
    }
    else{
        addLabels(glviewer);
        glviewer.render();
        modelLabels = true
    }
}

