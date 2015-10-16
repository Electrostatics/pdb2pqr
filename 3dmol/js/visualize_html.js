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
    "<br><font style='color:white; font-size:12pt'><center>Surface</font><hr width='75%'>" +
"    <select class='styled-select' id='selected_surface' onchange='update_surface(1)' style='max-width:50%;'>" +
"       <option style='color: black;' value='SAS'>Solvent Accessible</option>" +
"<option style='color: black;' value='SES'>Solvent Excluded</option>" +
"<option style='color: black;' value='VDW'>Van Der Waals</option>" +
"</select><br><br>" +

"<div class='inner'><ul class='button-group round' class='leftbutton'><input type='button' button class='button-backbone pure-button' style='width: 85px; height: 30px; color: black' value='On' onclick='on_surface()'></button></input>" +

"    <input type='button' button class='button-backbone pure-button' style='width: 85px; height: 30px; color: black' input type='button' value='Off' onclick='glviewer.removeSurface(surf)'></button></input></ul></div>" +

"<div class='inner'><ul class='button-group round' class='rightbutton'><input type='button' button class='button-labela pure-button' style='width: 85px; height: 30px; color: black; font-size: 10pt' value='Translucent' onclick='update_surface(2)'></button></input>" +

"    <input type='button' button class='button-unlabela pure-button' style='width: 85px; height: 30px; color: black' input type='button' value='Opaque' onclick='update_surface(3)'></button></input></ul></div>" +

"<font style='color:white; font-size:12pt'>Isosurface</font><hr width='75%'>" +

//change min isoval
" <p style='color:white; font-size: 16px'> Minimum kT/e <input type=range min=-50 max=50 value=-5 id='min_isoval2' step=1 onmouseup='set_min_isoval2(value)'>&nbsp;&nbsp;&nbsp;&nbsp;<span id='min_isoval'>-5 </span> </p>  " +

//change max isoval
" <p style='color:white; font-size: 16px'> Maximum kT/e <input type=range min=-50 max=50 value=5 id='max_isoval2' step=1 onmouseup='set_max_isoval2(value)'>&nbsp;&nbsp;&nbsp;&nbsp;<span id='max_isoval'>5 </span> </p>  " +

//reset button
"<div class='inner'><ul class='button-group round'></input>" +
"    <input type='button' button class='button-backbone pure-button' style='width: 85px; height: 30px; color: black' input type='button' value='Reset' onclick='reset_vals()'></button></input></ul></div><br><br>" +
"<font style='color:white; font-size:12pt'>Model</font><hr width='75%'>" +

"<select class='styled-select' id='selected_vis' onchange='set_vis()' style='max-width:50%;'>"+
"        <option style='color: black;' value='line'>Line </option> "+
"        <option style='color: black;' value='stick'>Stick </option>" +
"        <option style='color: black;' value='cross'>Cross </option>"+
"        <option style='color: black;' value='sphere'>Sphere </option>"+
"        <option style='color: black;' value='cartoon'>Cartoon </option>"+
"    </select>" 
+  "<br><br>"  

 //change color scheme 
 var c =   
"<font style='color:white; font-size:12pt'>Scheme </font><hr width='75%'><select class='styled-select' id='selected_scheme' onchange='update_surface(0)' style='max-width:50%;'>" +
"        <option style='color: black;' value='RWB'>Red-White-Blue </option>" +
"        <option style='color: black;' value='ROYGB'>Red-Green-Blue </option>" +
"        <!--<option style='color: black;' value='BWR'>Blue-White-Red </option>-->" +
"    </select><br><br>" 

var h =

"    <div class='inner'><ul class='button-group round'><input type='button' button class='button-labela pure-button' style='width: 85px; height: 30px; color: black' value='Add labels'" +
"        onclick='addLabels(glviewer); glviewer.render();'></button></input>" +
"    <input type='button' button class='button-unlabela pure-button' style='width: 85px; height: 30px; color: black' value='Remove labels'" +
"        onclick='removetheLabels(glviewer); glviewer.render();'></button></input></ul></div>"  +
"<div class='inner'><ul class='button-group round'><input type='button' button class='button-backbone pure-button' style='width: 85px; height: 30px; color: black' input type='button' value='Recenter' onclick='glviewer.zoomTo();'></button></input></ul></div>" +
"<br></div></div>" 





document.write(a)
document.write(b)
document.write(c)




document.write(h)

}