function build_page(jobid){
	//jobid = 1234;
	document.title = "3Dmol Visualization " + jobid

	var a = 
	"<div id='gldiv'></div>" +
	"<hr style='margin: 0;'>" +
	"<br>" +
    "<div id='outer'>" +
    "<div id='buttons'>" 
	
	//change surface
    var b = 
    "<font style='color: black; font-size: 8pt; font-style: italic'>Click to change</font><br>" +
"    <select class='styled-select' id='selected_surface' onchange='update_surface(1)' >" +
"       <option style='color: black;' value='SAS'>Solvent Accessible</option>" +
"<option style='color: black;' value='SES'>Solvent Excluded</option>" +
"<option style='color: black;' value='VDW'>Van Der Waals</option>" +
"</select>" +
"    <br><br>" 

 //change color scheme 
 var c =   
"<select class='styled-select' id='selected_scheme' onchange='update_surface(0)'>" +
"        <option style='color: black;' value='RWB'>Red-White-Blue </option>" +
"        <option style='color: black;' value='ROYGB'>Red-Green-Blue </option>" +
"        <option style='color: black;' value='BWR'>Blue-White-Red </option>" +
"    </select>" +
"<br> <br>"

//change min isoval
var d =
" <table><tr><td style='color: white; font-size: 18px'>Minimum: &nbsp;</td> " +
"<td><br><input type='text' id='min_isoval2' value='-5' style='text-align:right; width: 50px; height: 30px; background: transparent; border:0; color: white; font-size: 18px'> </td>" +
"<td style='color: white; font-size: 18px'>kT/e &nbsp;&nbsp;</td>" +
"<td valign='top'><br><br><input type=range min=-50 max=50 value=-5 id='min_isoval2' step=1 onmouseup='set_min_isoval2(value)'> </td></tr></table>"

//change max isoval
var e =
" <table><tr><td style='color: white; font-size: 18px'>Maximum: &nbsp;</td> " +
"<td><input type='text' id='max_isoval2' value='5' style='text-align:right; width: 50px; height: 30px; background: transparent; border:0; color: white; font-size: 18px'> </td>" +
"<td style='color: white; font-size: 18px'>kT/e &nbsp;&nbsp;</td>" +
"<td valign='top'><br><input type=range min=-50 max=50 value=5 id='max_isoval2' step=1 onmouseup='set_max_isoval2(value)'> </td></tr></table>"

//buttons
var f =
"<div class='inner'><ul class='button-group round'></input>" +
"    <input type='button' button class='button-backbone pure-button' style='width: 100px; height: 40px' input type='button' value='Reset' onclick='reset_vals()'></button></input></ul></div>" +
"    <br>" 

var g =
"  <div class='inner'><button id='visualbutton' id='leftbutton' data-dropdown='visualdrop' aria-controls='drop' aria-expanded='false' class='round button dropdown'>Visualization</button><br>" +
"<ul id='visualdrop' data-dropdown-content class='f-dropdown' role='menu' aria-hidden='false' tabindex='-1'>" +
"  <li><a onclick='glviewer.setStyle({},{stick:{}}); glviewer.render();'>Stick</a></li>" +
"  <li><a onclick='glviewer.setStyle({},{line:{}}); glviewer.render();'>Line</a></li>" +
"  <li><a onclick='glviewer.setStyle({},{cross:{linewidth:5}}); glviewer.render();'>Cross</a></li>" +
"  <li><a onclick='glviewer.setStyle({},{sphere:{}}); glviewer.render();'>Sphere</a></li>" +
"  <li><a onclick='glviewer.setStyle({hetflag:false},{cartoon:{color: 'spectrum'}}); glviewer.render();'>Cartoon</a></li>" +
"</ul></div>" +
"    <br>"

var h =
"<div class='inner'><ul class='button-group round' class='leftbutton'><input type='button' button class='button-backbone pure-button' style='width: 100px; height: 40px' value='On' onclick='on_surface()'></button></input>" +
"    <input type='button' button class='button-backbone pure-button' style='width: 100px; height: 40px' input type='button' value='Off' onclick='glviewer.removeSurface(surf)'></button></input></ul></div>" +
"    <div class='inner'><ul class='button-group round' class='rightbutton'><input type='button' button class='button-labela pure-button' style='width: 100px; height: 40px' value='Translucent' onclick='update_surface(2)'></button></input>" +
"    <input type='button' button class='button-unlabela pure-button' style='width: 100px; height: 40px' input type='button' value='Opaque' onclick='update_surface(3)'></button></input></ul></div>" +
"    <div class='inner'><ul class='button-group round'><input type='button' button class='button-labela pure-button' style='width: 100px; height: 40px' value='Add labels'" +
"        onclick='addLabels(glviewer); glviewer.render();'></button></input>" +
"    <input type='button' button class='button-unlabela pure-button' style='width: 100px; height: 40px' value='Remove labels'" +
"        onclick='removetheLabels(glviewer); glviewer.render();'></button></input></ul></div>" + "<div class='inner'><ul class='button-group round'><input type='button' button class='button-distance pure-button' style='width: 100px; height: 40px' value='' onclick='glviewer.zoomTo();'>" +
"    </button></input>" +
"   <input type='button' button class='button-recenter pure-button' style='width: 100px; height: 40px' value='Recenter' onclick='glviewer.zoomTo();'></button> </input></ul></div>" +
"<br></div></div>" 




document.write(a)
document.write(b)
document.write(c)
document.write(d)
document.write(e)
document.write(f)
document.write(g)
document.write(h)
}