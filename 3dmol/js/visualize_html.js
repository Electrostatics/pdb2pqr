function build_page(jobid){
	//jobid = 1234;
	document.title = "3Dmol Visualization " + jobid

	var a = 
	"<div id='gldiv'></div>" +
	"<!--<hr style='margin: 0;'>-->" +
	"<br>" +
    "<div id='outer'>" +
    "<div id='buttons'>" 
	
	//change surface
    var b = 
    "<!--<br><table border='1' width='100%' cellspacing='0' cellpadding='0'><tr><td valign='top'>--> " +
    "<font style='color:white; font-size:12pt'>Surface:</font>" + 
    "<br>" +
"    <select class='styled-select' id='selected_surface' onchange='update_surface(1)' style='max-width:50%;'>" +
"       <option style='color: black;' value='SAS'>Solvent Accessible</option>" +
"<option style='color: black;' value='SES'>Solvent Excluded</option>" +
"<option style='color: black;' value='VDW'>Van Der Waals</option>" +
"</select><br>" +

"<!--<div class='inner'><ul class='button-group round' class='leftbutton'><input type='button' button class='button-backbone pure-button' style='width: 85px; height: 30px; color: black' value='On' onclick='on_surface()'></button></input>" +

"    <input type='button' button class='button-backbone pure-button' style='width: 85px; height: 30px; color: black' input type='button' value='Off' onclick='glviewer.removeSurface(surf)'></button></input></ul></div><br>-->" +

"<table border='0' cellspacing='0' cellpadding='0'><tr><td>" +
"<label class='switch1'> " +
	"<input class='switch1-input' type='checkbox' onclick='surface_vis(this)' checked/>" +
	"<span class='switch1-label' data-on='Show' data-off='Hide'></span> " +
	"<span class='switch1-handle'></span>" +
"</label>" +

"</td>" +

"<td>" +
"<!--<div class='inner'><ul class='button-group round' class='rightbutton'><input type='button' button class='button-labela pure-button' style='width: 85px; height: 30px; color: black; font-size: 10pt' value='Translucent' onclick='update_surface(2)'></button></input>" +

"    <input type='button' button class='button-unlabela pure-button' style='width: 85px; height: 30px; color: black' input type='button' value='Opaque' onclick='update_surface(3)'></button></input></ul></div>-->" +

"<label class='switch2'>" +
	"<input class='switch2-input' type='checkbox' onclick='surface_opacity(this)' checked/>" +
	"<span class='switch2-label'   data-on='Opaque' data-off='Translucent'></span> " +
	"<span class='switch2-handle'></span> " +
"</label>" +

"</td></tr></table>"+

"<br><font style='color:white; font-size:12pt'>Isosurface:</font>" +

//change min isoval
" <p style='color:white; font-size: 16px'> Min<input type=range min=-50 max=50 value=-5 id='min_isoval2' step=1 oninput='set_min_isoval2(value)'>&nbsp;&nbsp;&nbsp;&nbsp;<span id='min_isoval'>-5 </span> kT/e </p>  " +

//change max isoval
" <p style='color:white; font-size: 16px'> Max<input type=range min=-50 max=50 value=5 id='max_isoval2' step=1 oninput='set_max_isoval2(value)'>&nbsp;&nbsp;&nbsp;&nbsp;<span id='max_isoval'> 5 </span> kT/e </p>  " +

//reset button
"<div class='inner'><ul class='button-group round'></input>" +
"    <input type='button' button class='button-backbone pure-button' style='width: 85px; height: 30px; color: black' input type='button' value='Reset' onclick='reset_vals()'></button></input></ul></div>" +

"<br><font style='color:white; font-size:12pt'>Model:</font>&nbsp;&nbsp;" +

"<select class='styled-select' id='selected_vis' onchange='set_vis()' style='max-width:50%;'>"+
"        <option style='color: black;' value='line'>Line </option> "+
"        <option style='color: black;' value='stick'>Stick </option>" +
"        <option style='color: black;' value='cross'>Cross </option>"+
"        <option style='color: black;' value='sphere'>Sphere </option>"+
"        <option style='color: black;' value='cartoon'>Cartoon </option>"+
"    </select>" 
+  "<br><br>"  +

 //change color scheme 
"<font style='color:white; font-size:12pt'>Scheme: </font>" +


 
"<select class='styled-select' id='selected_scheme' onchange='update_surface(0);show_colorbar();' style='max-width:50%;'>" +

"        <option style='color: black;' value='RWB'>Red-White-Blue </option>" +
"        <option style='color: black;' value='RGB'>Red-Green-Blue </option>" +
"        <!--<option style='color: black;' value='BWR'>Blue-White-Red </option>-->" +
"    </select><br>" +
"<!--<table border='0' cellspacing='0'><tr><td valign='top'>RWB<br>RGB</td><td valign='top'><img src='3dmol/images/rwb.png' width='250'><br><img src='3dmol/images/rgb.png' width='250'></td></tr></table>-->" +
 "<span id='colorbar'><img id='rwb' src='3dmol/images/rwb.png' width='250'></span>" +

"    <!--<div class='inner'><ul class='button-group round'><input type='button' button class='button-labela pure-button' style='width: 90px; height: 30px; color: black' value='Add labels'" +
"        onclick='addLabels(glviewer); glviewer.render();'></button></input>" +
"    <input type='button' button class='button-unlabela pure-button' style='width: 95px; height: 30px; color: black' value='Remove labels'" +
"        onclick='removetheLabels(glviewer); glviewer.render();'></button></input></ul></div>-->"  +

"<table><tr><td>" +
"<font style='color:white; font-size:12pt'>Labels: </font></td><td><label class='switch3'>" +
	"<input class='switch3-input' onclick='surface_labels(this)' type='checkbox' checked/>" +
	"<span class='switch3-label' data-on='Remove' data-off='Add'></span>" +
	"<span class='switch3-handle'></span> " +
"</label>" +

"</td></tr></table>"+


"<div class='inner'><ul class='button-group round'><input type='button' button class='button-backbone pure-button' style='width: 85px; height: 30px; color: black' input type='button' value='Recenter' onclick='glviewer.zoomTo();'></button></input></ul></div>" +
"<br></div></div>" +
"<!--</td></tr></table>-->" 





document.write(a)
document.write(b)
}