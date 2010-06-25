/*
LookSeq browser interface
(c) 2008 by Magnus Manske (mm6@sanger.ac.uk)
Released under GPL

    Petr Danecek (pd3@sanger.ac.uk), Team 145
*/

//_________________________________________________________________________________________________
// Global variables
//_________________________________________________________________________________________________
var cur_from = 0 ;
var cur_to = 0 ;
var is_loading = 0 ;
var chromosomes = new Array () ;
var chromosome_length = new Array () ;
var display_mode = 'indel' ;
var max_window = 1000;  // The maximum window - modified in runtime, value depends on the value of zoom set
var min_window = 110;   // Do not display less than this many bases - modified in runtime, value depends on image width which can be changed
var last_main_image_url = '' ;
var last_url_second = '' ;
var last_gc_url = '';
var last_coverage_url = '';
var last_annotation_url = '';
var lane_select = 'single' ;
var second_track_lanes = '' ;
var hbar_shown = 0 ;
var viewportwidth;
var viewportheight;
var search_results = new Array () ;
var fading = false ;
var ie = document.all;
var nn6 = document.getElementById && !document.all ;
var isdrag = false ;
var drag_x , drag_y ;
var dobj;

/*
Two functions not implemented in this script may be called. You can add them to a "custom.js" file in this directory if you wish.
function initialize_organism () {}
function single_position_double_click ( pos_x , pos_y ) {}
*/




// Hack around Internet Explorer 6.0 and earlier
if( typeof XMLHttpRequest == "undefined" ) XMLHttpRequest = function() {
  try { return new ActiveXObject("Msxml2.XMLHTTP.6.0") } catch(e) {}
  try { return new ActiveXObject("Msxml2.XMLHTTP.3.0") } catch(e) {}
  try { return new ActiveXObject("Msxml2.XMLHTTP") } catch(e) {}
  try { return new ActiveXObject("Microsoft.XMLHTTP") } catch(e) {}
  throw new Error( "This browser does not support XMLHttpRequest." )
};


//_________________________________________________________________________________________________
// Helper functions
//_________________________________________________________________________________________________

// Helper function used in on_search() and annotation_clicked()
function update_viewport_size () {
	if (typeof window.innerWidth != 'undefined') {
		viewportwidth = window.innerWidth,
		viewportheight = window.innerHeight
	} else if (typeof document.documentElement != 'undefined' && typeof document.documentElement.clientWidth != 'undefined' && document.documentElement.clientWidth != 0) {
		viewportwidth = document.documentElement.clientWidth,
		viewportheight = document.documentElement.clientHeight
	} else {
		viewportwidth = document.getElementsByTagName('body')[0].clientWidth,
		viewportheight = document.getElementsByTagName('body')[0].clientHeight
	}
}

// Used in fading
function grayOut(vis, options) {
  // Pass true to gray out screen, false to ungray
  // options are optional.  This is a JSON object with the following (optional) properties
  // opacity:0-100         // Lower number = less grayout higher = more of a blackout 
  // zindex: #             // HTML elements with a higher zindex appear on top of the gray out
  // bgcolor: (#xxxxxx)    // Standard RGB Hex color code
  // grayOut(true, {'zindex':'50', 'bgcolor':'#0000FF', 'opacity':'70'});
  // Because options is JSON opacity/zindex/bgcolor are all optional and can appear
  // in any order.  Pass only the properties you need to set.
  var options = options || {}; 
  var zindex = options.zindex || 50;
  var opacity = options.opacity || 70;
  var opaque = (opacity / 100);
  var bgcolor = options.bgcolor || '#000000';
  var dark=document.getElementById('darkenScreenObject');
  if (!dark) {
    // The dark layer doesn't exist, it's never been created.  So we'll
    // create it here and apply some basic styles.
    // If you are getting errors in IE see: http://support.microsoft.com/default.aspx/kb/927917
    var tbody = document.getElementsByTagName("body")[0];
    var tnode = document.createElement('div');           // Create the layer.
        tnode.style.position='absolute';                 // Position absolutely
        tnode.style.top='0px';                           // In the top
        tnode.style.left='0px';                          // Left corner of the page
        tnode.style.overflow='hidden';                   // Try to avoid making scroll bars            
        tnode.style.display='none';                      // Start out Hidden
        tnode.id='darkenScreenObject';                   // Name it so we can find it later
    tbody.appendChild(tnode);                            // Add it to the web page
    dark=document.getElementById('darkenScreenObject');  // Get the object.
  }
  if (vis) {
    // Calculate the page width and height 
    if( document.body && ( document.body.scrollWidth || document.body.scrollHeight ) ) {
        var pageWidth = document.body.scrollWidth+'px';
        var pageHeight = document.body.scrollHeight+'px';
    } else if( document.body.offsetWidth ) {
      var pageWidth = document.body.offsetWidth+'px';
      var pageHeight = document.body.offsetHeight+'px';
    } else {
       var pageWidth='100%';
       var pageHeight='100%';
    }   
	pageWidth='100%';
	pageHeight='100%' ;
    //set the shader to cover the entire page and make it visible.
    dark.style.opacity=opaque;                      
    dark.style.MozOpacity=opaque;                   
    dark.style.filter='alpha(opacity='+opacity+')'; 
    dark.style.zIndex=zindex;        
    dark.style.backgroundColor=bgcolor;  
    dark.style.width= pageWidth;
    dark.style.height= pageHeight;
    dark.style.display='block';				 
  } else {
     dark.style.display='none';
  }
}

function setOpacity(obj, opacity) {
  opacity = (opacity == 100)?99.999:opacity;

  if ( typeof obj == 'undefined' ) return ;
  if ( !obj ) return ;

  // IE/Win
  if ( is_MSIE() ) obj.style.filter="progid:DXImageTransform.Microsoft.Alpha(opacity=" + opacity + ");";

  // Safari<1.2, Konqueror
  obj.style.KHTMLOpacity = opacity/100;

  // Older Mozilla and Firefox
  obj.style.MozOpacity = opacity/100;

  // Safari 1.2, newer Firefox and Mozilla, CSS3
  obj.style.opacity = opacity/100;
}

function image_event_handler_get_position ( e , o ) {
  var ret = new Array () ;
  ret.x = 0 ;
  ret.y = 0 ;
	if (e.pageX || e.pageY) 	{
		ret.x = e.pageX;
		ret.y = e.pageY;
	}
	else if (e.clientX || e.clientY) 	{
		ret.x = e.clientX + document.body.scrollLeft + document.documentElement.scrollLeft;
		ret.y = e.clientY + document.body.scrollTop + document.documentElement.scrollTop;
	}
	
	try {
    ret.x -= o.offsetLeft + o.offsetParent.offsetLeft ;
    ret.y -= o.offsetTop + o.offsetParent.offsetTop ;
	} catch ( err ) {
    ret.x = -10000 ;
    ret.y = -10000 ;
	}
	
  return ret ;
}

// Hack around broken Internet Explorer positioning
function getElementPosition(elemID) {
    var offsetTrail = document.getElementById(elemID);
    var offsetLeft = 0;
    var offsetTop = 0;
    while (offsetTrail) {
        offsetLeft += offsetTrail.offsetLeft;
        offsetTop += offsetTrail.offsetTop;
        offsetTrail = offsetTrail.offsetParent;
    }
    if (navigator.userAgent.indexOf("Mac") != -1 && 
        typeof document.body.leftMargin != "undefined") {
        offsetLeft += document.body.leftMargin;
        offsetTop += document.body.topMargin;
    }
    return {left:offsetLeft, top:offsetTop};
}

function removeChildrenFromNode(node) {
   var len = node.childNodes.length;

	while (node.hasChildNodes())
	{
	  node.removeChild(node.firstChild);
	}
}

function is_MSIE () {
  return navigator.userAgent.indexOf ( "MSIE" ) > -1 ;
}


//_________________________________________________________________________________________________
// Here starts the fun
//_________________________________________________________________________________________________

function create_lane_items ( form , type ) {
    var tbl     = document.createElement("table");
    var tblBody = document.createElement("tbody");
	for ( var i = 0 ; i < lanes.length ; i++ ) {
		var name = new Array() ;
		var n ;
		var u = lanes[i].split('|') ;
		if ( u.length == 1 ) 
        {
            // Remove the .bam part of the name

			name = lanes[i];
            n = lanes[i].replace(/.bam$/,'').replace(/_/,'/');
		} else {
			name = u[0] ;
			name = name.match ( /^[^\.]+/ ) ;
			n = u[1] ;
		}
		var id = 'id4label_' + name ;
		var inp = document.createElement ( 'input' ) ;
		inp.id = id ;
		inp.type = type ;
		inp.name = 'chosenlane' ;
		inp.value = name ;
		inp.onclick = function () { activate_lane ( this.value ) ; } ;

		//  var label = document.createElement ( 'label' )
		//  if ( is_MSIE() ) label.setAttribute ( 'htmlFor' , id ) ;
		//  else label.setAttribute ( 'for' , id ) ;
		//  label.appendChild ( document.createTextNode ( n ) ) ;
		//  form.appendChild ( inp ) ;
		//  form.appendChild ( label ) ;
		//  form.appendChild ( document.createElement ( 'br' ) ) ;

        var row  = document.createElement("tr");
        var cell = document.createElement("td");
        cell.appendChild(inp);
        cell.appendChild(document.createTextNode(n));
        row.appendChild(cell);
        if ( typeof depth!='undefined' && depth[name] ) 
        {
            cell = document.createElement("td");
            var text = document.createTextNode('(' + depth[name] + 'X)');
            cell.appendChild(text);
            cell.style.fontSize = 'smaller';
            row.appendChild(cell);
        }
        tblBody.appendChild(row);
	}
    tbl.appendChild(tblBody);
    tbl.id = 'llist';
    form.appendChild(tbl);
}

function create_toggle_button ( form ) {
	// Toggle single/multi-lane mode
	var b = document.createElement ( 'input' ) ;
	b.id = 'button_single_multi' ;
	b.type = 'button' ;
	b.value = i18n['button_multiple_lanes'] ;
	b.onclick = toggle_lanes_radio_checkbox ;
	var msg = document.createElement ( 'div' ) ;
	msg.id = 'msg_update_manually' ;
	msg.style.display = 'none' ;
	msg.style.color = 'red' ;
	msg.style.fontSize = '7pt' ;
	msg.appendChild ( document.createTextNode ( i18n['daad1'] ) ) ;
	msg.appendChild ( document.createElement ( 'br' ) ) ;
	msg.appendChild ( document.createTextNode ( i18n['daad2'] ) ) ;
    // Uh, how not to add the button and stay compatible??
    b.style.display = 'none';
	form.appendChild ( b ) ;
	form.appendChild ( msg ) ;
	
	document.getElementById('lanes').appendChild ( form ) ;
}

function show_lanes () {
	var form = document.createElement ( 'form' ) ;
	var h = document.createElement ( 'p' ) ;
	form.name = 'lanes' ;
	form.id = 'lanes_form' ;
	h.id = 'lanes_title' ;
	h.style.fontWeight = 'bold' ;
	h.style.borderBottom = '1px solid black' ;
	h.appendChild ( document.createTextNode ( i18n['lanes_title'] ) ) ;
    // Do the title in lookseq.html instead.
	//  form.appendChild ( h ) ;
	var lanelist = document.createElement ( 'div' ) ;
	lanelist.id = "lanelist" ;
    // Do the style in .css instead.
    //
	//  lanelist.style.fontSize = "8pt" ;
	//  lanelist.style.maxHeight = "620px" ;
	//  lanelist.style.overflowY = 'scroll' ;
	//  lanelist.style.width = '100%' ;
	
	create_lane_items ( lanelist , 'radio' ) ;
	form.appendChild ( lanelist ) ;
	create_toggle_button ( form ) ;
}

// Event handler : toggles between single-lane and multi-lane mode, except for IE
function toggle_lanes_radio_checkbox () {
	for ( f = 0 ; f < document.forms.length ; f++ ) {
		if ( document.forms[f].name == 'lanes' ) break
	}
	var state = document.forms[f].elements[0].type ;
	var newstate = state == 'radio' ? 'checkbox' : 'radio' ;
	var firstid = '' ;
	lane_select = newstate == 'radio' ? 'single' : 'multiple' ;

	if ( is_MSIE() ) { // Did I mention I hate IE???
		var form = document.getElementById('lanes_form') ;
		var lt = document.getElementById('lanes_title') ;
		while ( lt.nextSibling ) {
			var i = lt.nextSibling ;
			if ( i.type == state && firstid == '' && i.checked ) {
				firstid = i.id ;
			}
			form.removeChild ( i ) ;
		}
		create_lane_items ( form , newstate ) ;
		create_toggle_button ( form ) ;
	} else {
		for ( var n = 0 ; n < document.forms[f].elements.length ; n++ ) {
			if ( document.forms[f].elements[n].type != state ) continue ;
			if ( firstid == '' && document.forms[f].elements[n].checked ) firstid = document.forms[f].elements[n].id ;
			document.forms[f].elements[n].type = newstate ;
		}
	}

	var b = document.getElementById ( 'button_single_multi' ) ;
	b.value = newstate == 'radio' ? i18n['button_multiple_lanes'] : i18n['button_single_lane'] ;
	document.getElementById('msg_update_manually').style.display = 'none' ;
	if ( firstid ) {
		document.getElementById(firstid).checked = true ;
		document.getElementById(firstid).onclick () ;
		update_image();
	}
}

// Event handler : User selected a new lane
function activate_lane ( new_lane ) {
	if ( new_lane == current_lane && lane_select == 'single' ) return ;
	if ( lane_select == 'multiple' ) { // No live update when using multiple lanes
		document.getElementById('msg_update_manually').style.display = 'block' ;
		return ;
	}
	old_lane = current_lane ;
	current_lane = new_lane ;

	if ( is_MSIE() ) {
		if ( old_lane ) document.getElementById('id4label_'+old_lane).checked = false ;
		document.getElementById('id4label_'+current_lane).checked = true ;
	}
	
	update_image () ;
}

// Updates the reference link. Called from update_image()
function update_reflink ( lanes , display ) {
	var sel = document.getElementById("chr_list");
	var reflink = cgi_path + '/index.pl?show=' ;
	reflink += sel.options[sel.selectedIndex].value + ':' + cur_from + '-' + cur_to + ',' + display_mode ;
	reflink += '&lane=' + lanes + '&width=' + img_width ;
    reflink += '&win=' + max_window;
	reflink += '&display=' + display ;
	if ( indel_zoom != 'auto' ) reflink += '&maxdist=' + indel_zoom ;
    if ( mapq_cutoff != 0 ) reflink += '&mapq=' + mapq_cutoff;
	if ( document.getElementById('display_second_track').checked ) {
		reflink += '&second_image=' + second_track_lanes ;
		if ( document.getElementById('squeeze_tracks').checked )
			reflink += '&squeeze_tracks=1' ;
	}
	document.getElementById('reflink').href = reflink ;
}

function get_max_to(cur_from)
{
    var max_to = cur_from + max_window;
    var chr = get_selected_chromosome();
    if ( max_to > chromosome_length[chr] )  max_to = chromosome_length[chr];
    return max_to;
}

// Updates main image and annotation/GC display, according to current settings. The heart of this script.
function update_image () 
{
	if ( is_loading > 0 ) 
    {
        // The user pressed the refresh button while 'Updating..' in progress. This usually
        //  happens when get_data.pl or other script dies. In such a case, the script was dead
        //  and the user was left with the only option: reload of the page. The new behaviour:
        //  If the page is not responding 15 sec, the user can change the region.

        var end_time = new Date().getTime();
        if ( end_time-start_time<15*1000 )
        {
            setTimeout(update_image , 100);
            return ;
        }
        loading(-is_loading);
	}
    start_time = new Date().getTime();

    // This seems to be "nice" spacing of chars for the given image width.
    document.getElementById("image_width").value = img_width;
    min_window = Math.floor(110*(img_width/700.));

    // In our application, we cannot allow a too large window size. Let the size be fixed,
    //  modifiable only by the buttons.
	cur_from   = Math.floor ( document.getElementById('chr_from').value.replace(/,/g,'') ) ;
    var max_to = chromosome_length[get_selected_chromosome()];
    cur_to = get_max_to(cur_from);
    if ( cur_from > cur_to )
    {
        cur_from = cur_to - min_window;
        document.getElementById('chr_from').value = cur_from;
    }
    document.getElementById('chr_to').value = cur_to ;

    if ( display_mode!='indel' && cur_to-cur_from>min_window*1.1 )
        cur_to = cur_from+min_window;

	//  if ( 0 && cur_to - cur_from > min_window*1.1 ) {
	//  	if ( display_mode != 'indel' ) {
	//  		document.getElementById("display_mode_indel").checked = true ;
	//  		document.getElementById("display_mode_indel").selected = true ;
	//  		display_mode = 'indel' ;
    //          window.alert("switching to indel mode: "+cur_to+" "+cur_from+" "+min_window);
	//  	}
	//  	document.getElementById("display_mode_pileup").disabled = true ;
	//  	document.getElementById("display_mode_paired_pileup").disabled = true ;
	//  } else {
	//  	document.getElementById("display_mode_pileup").disabled = false ;
	//  	document.getElementById("display_mode_paired_pileup").disabled = false ;
	//  }

    if ( display_mode == 'paired_pileup' )
        document.getElementById("display_mode_paired_pileup").selected = true ;
    if ( display_mode == 'indel' )
        document.getElementById("display_mode_indel").selected = true ;

	document.getElementById("legend_indel").style.display = display_mode == 'indel' ? 'block' : 'none' ;
	document.getElementById("legend_pileup").style.display = display_mode == 'pileup' || display_mode == 'paired_pileup' ? 'block' : 'none' ;
	document.getElementById("legend_coverage").style.display = display_mode == 'coverage' ? 'block' : 'none' ;
	document.getElementById("legend_known_snps").style.display = document.getElementById('display_known_snps').checked ? 'block' : 'none' ;
	document.getElementById("legend_annotation").style.display = document.getElementById('display_annotation').checked ? 'block' : 'none' ;
	document.getElementById("legend_gc").style.display = document.getElementById('display_gc').checked ? 'block' : 'none' ;
	document.getElementById("legend_inv").style.display = document.getElementById('display_inversions_ext').checked ? 'block' : 'none' ;
	document.getElementById("legend_coverage_two_samples").style.display = document.getElementById('display_coverage').checked && document.getElementById('display_second_track').checked ? 'block' : 'none' ;
	
	
	document.getElementById("display_pair_links").disabled = display_mode != 'indel' && display_mode != 'paired_pileup' ;
	document.getElementById('display_inversions_ext').disabled = document.getElementById('display_inversions').checked && display_mode == 'indel' ? false : true ;
	document.getElementById('msg_update_manually').style.display = 'none' ;

	//  document.getElementById("button_zoom_150").disabled = cur_to - cur_from + 1 <= min_window+2 ;
	//  document.getElementById("button_zoom_in").disabled = cur_to - cur_from + 1 <= min_window+2 ;
	//  document.getElementById("button_zoom_1500").disabled = cur_to - cur_from == 2000 ;
	//  document.getElementById("button_zoom_50k").disabled = cur_to - cur_from == 50000 ;
	//  document.getElementById('button_left').disabled = cur_from == 1 ;
	//  document.getElementById('button_right').disabled = cur_to == max_to ;
	//  document.getElementById('button_zoom_out').disabled = cur_from == 1 && cur_to == max_to ;
	//  document.getElementById('button_zoom_chr').disabled = cur_from == 1 && cur_to == max_to ;

	var show_annotation = document.getElementById('display_annotation').checked ? 1 : 0 ;
	var show_gc = document.getElementById('display_gc').checked ? 1 : 0 ;
	var show_coverage = document.getElementById('display_coverage').checked ? 1 : 0 ;
	
	var display = '|' ;
	display += document.getElementById('display_perfect').checked ? 'perfect|' : '' ;
	display += document.getElementById('display_snps').checked ? 'snps|' : '' ;
	display += document.getElementById('display_single_reads').checked ? 'single|' : '' ;
	display += document.getElementById('display_inversions').checked ? 'inversions|' : '' ;
	display += document.getElementById('display_inversions_ext').checked ? 'inversions_ext|' : '' ;
	display += document.getElementById('display_pair_links').checked ? 'pairlinks|' : '' ;
	display += document.getElementById('display_known_snps').checked ? 'potsnps|' : '' ;
	display += document.getElementById('display_uniqueness').checked ? 'uniqueness|' : '' ;
	display += document.getElementById('display_gc').checked ? 'gc|' : '' ;
	display += document.getElementById('display_coverage').checked ? 'coverage|' : '' ;
	display += document.getElementById('show_arrows').checked ? 'orientation|' : '';
	
	document.getElementById('annotation_image').style.display = show_annotation ? 'block' : 'none' ;
	document.getElementById('gc_image').style.display = show_gc ? 'block' : 'none' ;
	document.getElementById('coverage_image').style.display = show_coverage ? 'block' : 'none' ;
	
	var sel = document.getElementById("chr_list");
    indel_zoom = document.getElementById('indel_zoom').value ;
    mapq_cutoff = document.getElementById('mapq') ? document.getElementById('mapq').value : 0;

	var url = cgi_path + '/get_data.pl' ;
	url += '?from=' + cur_from ;
	url += '&to=' + cur_to ;
	url += '&chr=' + sel.options[sel.selectedIndex].value ;
	url += '&output=image' ;
	url += '&width=' + img_width ;
	var lanes = condense_lanes_for_url () ;
	url += '&lane=' + lanes ;
    //if ( !document.getElementById('use_cache').checked ) 
    //    url += '&clean=1';
	var annotation_url = url + '&view=annotation' ;
	var gc_url = url + '&view=gc' ;

	var url_second = url ;
	var coverage_url = url + '&view=coverage&display=|noscale' + display ;

	var url_part2 = '' ;
	url_part2 += '&view=' + display_mode ;
	if ( indel_zoom != 'auto' ) url_part2 += '&maxdist=' + indel_zoom ;
    if ( mapq_cutoff != 0 ) url_part2 += '&mapq=' + mapq_cutoff;
	url_part2 += '&display=' + display ;
	url_part2 += '&debug=' ;
	url_part2 += test ? '1' : '0' ;
	if ( document.getElementById('squeeze_tracks').checked ) url_part2 += '&height=256' ;
	url += url_part2 ;
	
	document.getElementById('vruler').style.display = 'none' ;
	
	document.getElementById("link_to_text").style.display = display_mode == 'pileup' ? 'inline' : 'none' ;
	if ( display_mode == 'pileup' || display_mode == 'paired_pileup' ) document.getElementById("link_to_text").href = url.replace('&output=image','&output=text') ;
	
	display += show_annotation ? 'annotation|' : '' ; // Is ignored by image drawing routine, but used by update_reflink()
	display += show_gc ? 'gc|' : '' ; // Is ignored by image drawing routine, but used by update_reflink()
	display += show_coverage ? 'coverage|' : '' ; // Is ignored by image drawing routine, but used by update_reflink()
	update_reflink ( lanes , display ) ;

	var dst = document.getElementById('display_second_track').checked ;
	var load = 0 ;
	var opaq = 20 ;
	if ( last_main_image_url != url ) {
		setOpacity ( document.getElementById('main_image') , opaq ) ;
		document.getElementById('main_image').src = url ;
		load++ ;
		last_main_image_url = url ;
	}
	if ( dst && last_url_second != url_second ) {
		url_second += '&lane=' + second_track_lanes + url_part2 ;
		setOpacity ( document.getElementById('second_image') , opaq ) ;
		document.getElementById('second_image').src = url_second ;
		load++ ;
		last_url_second = url_second ;
	}
	if ( show_annotation && last_annotation_url!=annotation_url ) 
    {
		setOpacity ( document.getElementById('annotation_image') , opaq ) ;
		document.getElementById('annotation_image').src = annotation_url ;
		load++ ;
        last_annotation_url = annotation_url;
	}
	if ( show_gc && last_gc_url!=gc_url ) 
    {
        // If the URL does not change, Safari browser will never fire onload event.
		setOpacity ( document.getElementById('gc_image') , opaq ) ;
		document.getElementById('gc_image').src = gc_url ;
		load++ ;
        last_gc_url = gc_url;
	}
	if ( show_coverage && last_coverage_url!=coverage_url) 
    {
		setOpacity ( document.getElementById('coverage_image') , opaq ) ;
		if ( dst ) coverage_url += "&second=" + second_track_lanes ;
		document.getElementById('coverage_image').src = coverage_url ;
		load++ ;
        last_coverage_url = coverage_url;
	}
	update_hbar () ;
	loading ( load ) ;
}

function condense_lanes_for_url () {
	var lanes ;
	if ( lane_select == 'single' ) {
		lanes = current_lane ;
	} else {
		lanes = new Array () ;
		for ( f = 0 ; f < document.forms.length ; f++ ) {
			if ( document.forms[f].name == 'lanes' ) break
		}
		for ( n = 0 ; n < document.forms[f].elements.length ; n++ ) {
			if ( document.forms[f].elements[n].type != 'checkbox' ) continue ;
			if ( !document.forms[f].elements[n].checked ) continue ;
			lanes.push ( document.forms[f].elements[n].value ) ;
		}
		lanes = lanes.join ( ',' ) ;
	}
	return lanes ;
}


// Event handler : Button to move image to the left was clicked.
function image_left () {
	var diff = cur_to - cur_from + 1 ;
	diff = Math.floor ( diff / 4 ) ;
	if ( cur_from - diff < 1 ) {
		diff = cur_from - 1 ;
	}
	document.getElementById('chr_from').value = Math.floor ( cur_from ) - Math.floor ( diff ) ;
	document.getElementById('chr_to').value = Math.floor ( cur_to ) - Math.floor ( diff ) ;
	update_image () ;
}

// Event handler : Button to move image to the right was clicked.
function image_right () {
	var diff = cur_to - cur_from + 1 ;
	diff = Math.floor ( diff / 4 ) ;
	document.getElementById('chr_from').value = Math.floor ( cur_from ) + Math.floor ( diff ) ;
	document.getElementById('chr_to').value = Math.floor ( cur_to ) + Math.floor ( diff ) ;
	update_image () ;
}

// After loading an image, set its properties back to normal
function reset_image_display ( id ) {
	var elm = document.getElementById ( id ) ;
	elm.style.left = 0 ;
	elm.style.width = '' ;
	elm.style.height = '' ;
	setOpacity ( elm , 100 ) ;
}

// Callback function for loading the main image.
function main_image_loaded () {
	//document.getElementById('img_container').style.width = img_width + 'px' ;
	//document.getElementById('legend_container').style.width = img_width + 'px' ;
	reset_image_display ( 'main_image' ) ;
	loading ( -1 ) ;
}

// Callback function for loading the secondary image.
function second_image_loaded () {
	//document.getElementById('img_container').style.width = img_width + 'px' ;
	//document.getElementById('legend_container').style.width = img_width + 'px' ;
	reset_image_display ( 'second_image' ) ;
	loading ( -1 ) ;
}

// Callback function for loading the annotation image.
function annotation_image_loaded () {
	reset_image_display ( 'annotation_image' ) ;
	loading ( -1 ) ;
}

// Callback function for loading the GC% image.
function gc_image_loaded () {
	reset_image_display ( 'gc_image' ) ;
	loading ( -1 ) ;
}

// Callback function for loading the coverage image.
function coverage_image_loaded () {
	reset_image_display ( 'coverage_image' ) ;
	loading ( -1 ) ;
}

// Start/end loading of images. Handles transparency issues.
function loading ( diff ) {
    //    window.alert('loading='+is_loading+' diff='+diff);
	if ( is_loading == 0 && diff > 0 ) { // Turn on loading
		document.getElementById('loading').style.display = 'inline' ;
	}
	is_loading += diff ;

    if ( is_loading<0 )
        is_loading = 0;

	if ( is_loading == 0 ) { // Turn off loading
		document.getElementById('loading').style.display = 'none' ;
		document.getElementById('img_container').style.position = '' ;
	}
}

// Event handler : Squeeze tracks toggled
function toggle_squeeze_tracks () {
	update_image () ;
}

// Event handler : Second track toggled
function toggle_second_track () {
	var checked = document.getElementById('display_second_track').checked ;
	if ( checked ) {
		second_track_lanes = condense_lanes_for_url () ;
		show_second_image () ;
		document.getElementById('squeeze_tracks').checked = true ;
		update_image () ;
	} else {
		document.getElementById('squeeze_tracks').checked = false ;
		document.getElementById('squeeze_tracks').disabled = true ;
		document.getElementById('second_image').style.display = 'none' ;
		document.getElementById('second_image_name').style.display = 'none' ;
		update_image () ;
	}
}

// Event handler : Main image was double-clicked
function main_image_clicked ( event ) {
	document.getElementById('img_container').style.position = '' ;

    var pos_x = event.offsetX?(event.offsetX):event.pageX-document.getElementById("main_image").offsetLeft;
    var pos_y = event.offsetY?(event.offsetY):event.pageY-document.getElementById("main_image").offsetTop;
	
	var new_middle = Math.floor ( Math.floor ( cur_to - cur_from + 1 ) * pos_x / document.getElementById("main_image").width ) + cur_from ;

	//var new_range = Math.floor ( ( cur_to - cur_from + 1 ) / 2 ) ;
	//if ( new_range < 100 ) 
    //{
	//	if ( window.single_position_double_click ) 
    //    {   
    //        single_position_double_click ( pos_x , pos_y ) ;
    //    }
    //    return;
	//}
	// show_new_range ( new_middle , new_range ) ;

	show_new_range ( new_middle , cur_to-cur_from+1) ;
}

// Switches the display to a new range.
function show_new_range ( new_middle , new_range ) 
{
    new_middle = parseInt(new_middle);
    new_range  = parseInt(new_range);

    // Not sure if this was a good thing to do
	//  if ( new_range + 1 > new_middle * 2 ) new_range = new_middle * 2 - 1 ;

    // Do not zoom more than this much
	if ( new_range < min_window ) new_range = min_window;
	
	var new_from = Math.floor ( new_middle - new_range / 2 ) ;
	var new_to = Math.floor ( new_from + new_range ) ;
	if ( new_from < 1 ) new_from = 1 ;
	if ( new_to < new_from + min_window ) new_to = new_from + min_window;
	
    max_window = new_to - new_from;

    // if ( display_mode!='indel' && new_to - new_from > min_window*1.1 ) 
    // {
    //     window.alert("would switch to indel mode:"+
    //         " new_from="+new_from+
    //         " new_to="+new_to+
    //         " min_win="+min_window+
    //         " max_win="+max_window+
    //         "");
    // }

	zoom_image ( new_from , new_to ) ;
	document.getElementById('chr_from').value = new_from ;
	document.getElementById('chr_to').value = new_to ;
	update_image () ;
}

// Event handler : Button to zoom out was clicked.
function zoom_out () {
	var new_middle = Math.floor ( cur_from + ( cur_to - cur_from ) / 2 ) ;
	var new_range = Math.floor ( ( cur_to - cur_from + 1 ) * 2 ) ;
	var new_from = Math.floor ( new_middle - new_range / 2 ) ;
	var new_to = new_from + new_range ;
	while ( new_from < 1 ) {
		new_from++ ;
		new_to++ ;
	}
	
	zoom_image ( new_from , new_to ) ;
	document.getElementById('chr_from').value = new_from ;
	document.getElementById('chr_to').value = new_to ;
	update_image () ;
}

// Stretches/comresses the image prior to loading the new view. For effect only...
function zoom_image ( new_from , new_to ) {
    // In my browser it does terrible things, better to disable this effect :)
    return;

	var dfrom = cur_from - new_from ;
	var dto = cur_to - new_to ;
	var ow = cur_to - cur_from + 1 ;
	var nw = new_to - new_from + 1 ;
	
	zoom_single_image ( 'main_image' , img_width * dfrom / nw , img_width * ow / nw ) ;
	if ( document.getElementById('display_gc').checked ) zoom_single_image ( 'gc_image' , img_width * dfrom / nw , img_width * ow / nw ) ;
	if ( document.getElementById('display_annotation').checked ) zoom_single_image ( 'annotation_image' , img_width * dfrom / nw , img_width * ow / nw ) ;
	if ( document.getElementById('display_coverage').checked ) zoom_single_image ( 'coverage_image' , img_width * dfrom / nw , img_width * ow / nw ) ;
}

function zoom_single_image ( id , left , width ) {
	var elm = document.getElementById ( id ) ;
	var h = elm.height ;
	elm.style.left = left + "px" ;
	elm.style.width = width + "px" ;
	elm.style.height = h + "px" ;
}

// Event handler : Button to zoom in was clicked.
function zoom_in () {
	var new_middle = Math.floor ( cur_from + ( cur_to - cur_from ) / 2 ) ;
	var new_range = Math.floor ( ( cur_to - cur_from + 1 ) / 2 ) ;
	if ( new_range < min_window ) new_range = min_window;
	var new_from = Math.floor ( new_middle - new_range / 2 ) ;
	var new_to = new_from + new_range ;
	while ( new_from < 1 ) {
		new_from++ ;
		new_to++ ;
	}
	
	zoom_image ( new_from , new_to ) ;
	document.getElementById('chr_from').value = new_from ;
	document.getElementById('chr_to').value = new_to ;
	update_image () ;
}

// Event handler : Mouse is moving over main image. display vertical marker bar if appropriate.
function mouse_moved_over_main_image ( event ) {
	if ( isdrag ) return ;
	
    var pos_x = event.offsetX?(event.offsetX):event.pageX - document.getElementById("main_image").offsetLeft;
    var pos_y = event.offsetY?(event.offsetY):event.pageY - document.getElementById("main_image").offsetTop;
	
	var pos = Math.floor ( Math.floor ( cur_to - cur_from + 1 ) * pos_x / document.getElementById("main_image").width ) + cur_from ;
	var text = document.createTextNode ( pos ) ;
	var md = document.getElementById ( 'mini_data' ) ;
	removeChildrenFromNode ( md ) ;
	md.appendChild ( text ) ;
	md.style.display = 'inline' ;
	
	pos_x += getElementPosition('main_image').left - 1 ;
	
	if ( is_loading > 0 ) return ;
	
	var h = document.getElementById("main_image").height ;
	h += document.getElementById('display_second_track').checked ? 2+document.getElementById('second_image').clientHeight : 0 ;
	h += document.getElementById('display_annotation').checked ? 2+document.getElementById('annotation_image').clientHeight : 0 ;
	h += document.getElementById('display_gc').checked ? 2+document.getElementById('gc_image').clientHeight : 0 ;
	h += document.getElementById('display_coverage').checked ? 2+document.getElementById('coverage_image').clientHeight : 0 ;

    if ( document.getElementById('img_block') )
        h = document.getElementById('img_block').clientHeight;

	var vruler = document.getElementById('vruler') ;
	vruler.style.height = h + "px" ;
	vruler.style.left = pos_x + "px" ;
	vruler.style.top = getElementPosition('main_image').top + "px" ;
	vruler.style.display = 'block' ;
}

// Hides the vertical marker.
function hide_vruler () {
	document.getElementById('vruler').style.display = 'none' ;
	document.getElementById('mini_data').style.display = 'none' ;
}

// Gets a list of chromosomes in the organism. Needs to be done only once.
function initalize_chromosomes () {
	var url = cgi_path + '/chromosome_stats.pl' ;
	var request =  new XMLHttpRequest();
    //window.alert("testing in progress, please return later .. url="+url);
	request.open("GET", url, false);
	request.send(null);
	var lines = request.responseText.split ( "\n" ) ;
	for ( var i = 0 ; i < lines.length ; i++ ) {
		var l = lines[i].split ( "\t" ) ;
		if ( l[0] == '' ) break ;
        if ( l[0] == 'max_window' ) 
        {
            max_window = Math.floor (l[1]);
            continue;
        }
		chromosomes.push ( l[0] ) ;
		chromosome_length[l[0]] = Math.floor ( l[1] ) ;
		var o = document.createElement ( 'option' ) ;
		o.value = l[0] ;
		o.appendChild ( document.createTextNode ( l[0] ) ) ;
		document.getElementById('chr_list').appendChild ( o ) ;
	}
}

// Returns currently selected chromosome name.
function get_selected_chromosome () {
	var sel = document.getElementById("chr_list");
	return sel.options[sel.selectedIndex].value ;
}

// Event handler : Resizes display to show entire chromosome.
function full_chromosome_rezise () {
	document.getElementById('chr_from').value = 1 ;
	document.getElementById('chr_to').value = get_max_to(1);
}

// Event handler : Chanages the chromosome, full display.
function chromosome_changed () {
	full_chromosome_rezise () ;
	update_image () ;
}

function show_arrows_changed()
{
    update_image();
}

function switch_display_mode(mode)
{
    if ( mode=='indel' ) 
        display_mode_indel();
    else if ( mode=='paired_pileup' )
        display_mode_paired_pileup();
    else
        window.alert("FIXME");
}

function set_mode_indel () {
	display_mode = 'indel' ;
	document.getElementById('indel_zoom').disabled = false ;
}


// Event handler : Changes to InDel display mode.
function display_mode_indel () {
    set_mode_indel();
	update_image () ;
}

// Event handler : Changes to pileup display mode.
function display_mode_pileup () {
	display_mode = 'pileup' ;
	document.getElementById('indel_zoom').disabled = true ;
	update_image () ;
}

function set_mode_paired_pileup()
{
    display_mode = 'paired_pileup' ;
	document.getElementById('indel_zoom').disabled = true ;
}

// Event handler : Changes to paired pileup display mode.
function display_mode_paired_pileup () {
    set_mode_paired_pileup();
	update_image () ;
}

// Event handler : Changes to cooverage display mode.
function display_mode_coverage () {
	display_mode = 'coverage' ;
	document.getElementById('indel_zoom').disabled = true ;
	update_image () ;
}

// Event handler : Zooms to 1:1 display.
function zoom_150 () {
	var new_middle = Math.floor ( ( cur_from + cur_to ) / 2 ) ;
	var new_range = min_window;
    set_mode_paired_pileup();
	show_new_range ( new_middle , new_range ) ;
}

// Event handler : Zooms to 50kbp display.
function zoom_50k () {
	var new_middle = Math.floor ( ( cur_from + cur_to ) / 2 ) ;
	var new_range = 50000 ;
	show_new_range ( new_middle , new_range ) ;
    display_mode_indel();
}

function zoom_100k () {
	var new_middle = Math.floor ( ( cur_from + cur_to ) / 2 ) ;
	var new_range = 100000 ;
	show_new_range ( new_middle , new_range ) ;
    display_mode_indel();
}

function zoom (num) {
	var new_middle = Math.floor ( ( cur_from + cur_to ) / 2 ) ;
	max_window = parseInt(num);
    if ( max_window < min_window )
        max_window = min_window;

    if ( max_window == min_window )
        set_mode_paired_pileup();
    else
        set_mode_indel();
	show_new_range ( new_middle , max_window ) ;
}

// Event handler : Zooms to 2kb display.
function zoom_1500 () {
	var new_middle = Math.floor ( ( cur_from + cur_to ) / 2 ) ;
	var new_range = 2000 ;
	show_new_range ( new_middle , new_range ) ;
    display_mode_indel();
}

// Event handler : Handles InDel view "depth" changes.
function indel_zoom_changed () {
	var new_indel_zoom = document.getElementById('indel_zoom').value ;
	if ( indel_zoom == new_indel_zoom ) return ;
	indel_zoom = new_indel_zoom ;
	update_image () ;
}

function mapq_cutoff_changed () {
	var new_mapq_cutoff = document.getElementById('mapq').value ;
	if ( mapq_cutoff == new_mapq_cutoff ) return ;
	mapq_cutoff = new_mapq_cutoff;
	update_image () ;
}

// Event handler : Image width has been changed.
function image_width_changed () {
	var old_img_width = img_width ;
	img_width = document.getElementById('image_width').value ;
    max_window = Math.floor(img_width*110./700);

	var new_middle = Math.floor ( ( cur_from + cur_to ) / 2 ) ;
	var new_range = Math.floor ( ( cur_to - cur_from ) * img_width / old_img_width ) ;
	show_new_range ( new_middle , new_range ) ;
}

// Initializes the display settings from URL.
function initialize_display () {
	if ( display_init == '' ) return ;
	var a = display_init.split ( '|' ) ;
	var b = new Array () ;
	for ( var i = 0 ; i < a.length ; i++ ) {
		b[a[i]] = 1 ;
	}
	document.getElementById('display_perfect').checked		= b['perfect'] ? true : false ;
	document.getElementById('display_snps').checked 		= b['snps'] ? true : false ;
	document.getElementById('display_single_reads').checked	= b['single'] ? true : false ;
	document.getElementById('display_inversions').checked	= b['inversions'] ? true : false ;
	document.getElementById('display_inversions_ext').checked	= b['inversions_ext'] ? true : false ;
	document.getElementById('display_pair_links').checked	= b['pairlinks'] ? true : false ;
	document.getElementById('display_known_snps').checked	= b['potsnps'] ? true : false ;
	document.getElementById('display_uniqueness').checked	= b['uniqueness'] ? true : false ;
	document.getElementById('display_annotation').checked	= b['annotation'] ? true : false ;
	document.getElementById('display_gc').checked			= b['gc'] ? true : false ;
	document.getElementById('display_coverage').checked		= b['coverage'] ? true : false ;
	
	if ( indel_zoom != 'auto' ) document.getElementById('indel_zoom').value = indel_zoom ;
}

// Sanger-specific function - not used, ignore
function toggle_sanger_stuff () {
	var new_state = document.getElementById('collapseHEAD').style.display == '' ? 'none' : '' ;
	document.getElementById('collapseHEAD').style.display = new_state ;
	document.getElementById('nav_bar').style.display = new_state ;
	document.getElementById('navblock').style.display = new_state ;
	toggle();
	document.getElementById('button_sanger_bars').value = new_state == '' ? "Hide Sanger bars" : "Show Sanger bars" ;
}

// Updates the bar that shows the position of the current region within the chromosome.
function update_hbar () {
	var hbar = document.getElementById ( 'hbar' ) ;
	var dragbar = document.getElementById ( 'dragbar' ) ;
	if ( !hbar_shown ) {
		hbar.style.display = 'block' ;
		hbar.style.width = img_width + 'px' ;
		hbar.style.height = '16px' ;
		dragbar.style.height = ( parseInt ( hbar.style.height ) - 4 ) + "px" ;
		hbar_shown = 1 ;
	}
	
	var curw   = cur_to - cur_from + 1 ;
    var max_to = get_max_to(0);
	dragbar.style.left = Math.floor ( img_width * cur_from / max_to ) + "px" ;
	dragbar.style.width = Math.floor ( img_width * curw / max_to ) + "px" ;
}

// Gradual fade in/out for feature searching/inspecting. Recursive calling via timeout.
function fade ( id , o , step ) {
	fading = true ;
	setOpacity ( document.getElementById ( id ) , o ) ;
	if ( step == -1 && o <= 0 ) {
		grayOut ( 0 ) ;
		fading = false ;
		return ;
	}
	if ( step == 1 && o >= 80 ) {
		fading = false ;
		return ;
	}
	setTimeout ( "fade('"+id+"',"+(o+step*10)+","+step+");" , 5 ) ;
}

// Event handler : Search button has been clicked.
function on_search () {
	grayOut ( 1 , {'opacity':'0'} ) ;
	setTimeout ( "fade('darkenScreenObject',0,1);" , 1 ) ;
	var s = document.getElementById('search_query').value ;
	var url = cgi_path + '/query.pl?search=' + escape ( s ) ;
	var request =  new XMLHttpRequest();
	request.open("GET", url, false);
	request.send(null);
	search_results = request.responseText.split ( "\n" ) ;
	
	if ( search_results[0] == 'ERROR' ) {
		if ( fading && !is_MSIE() ) fade ( 'darkenScreenObject' , 80 , 1 ) ;
		alert ( search_results[1] ) ;
		if ( is_MSIE() ) grayOut ( 0 ) ;
		else setTimeout ( "fade('darkenScreenObject',80,-1);" , 1 ) ;
		return ;
	}


	update_viewport_size () ;
	var darken = document.getElementById('darkenScreenObject') ;
	var centerx = Math.floor ( viewportwidth / 2 ) ;
	var centery = Math.floor ( viewportheight / 2 ) ;
	var w2 = 250 ; // Half the width of the white box
	var h2 = 200 ; // Half the height of the white box
	var d = document.createElement ( 'div' ) ;
	d.id = 'white_wizard' ;
	d.style.zindex = 51 ;
	d.style.position = 'absolute' ;
	d.style.left = ( centerx - w2 ) + "px" ;
	d.style.top = ( centery - h2 ) + "px" ;
	d.style.width = ( w2 * 2 ) + "px" ;
	d.style.height = ( h2 * 2 ) + "px" ;
	d.style.backgroundColor = 'white' ;
	
	var subd = document.createElement ( 'div' ) ;
	subd.style.width = ( w2 * 2 ) + "px" ;
	subd.style.height = ( h2 * 2 - 30 ) + "px" ;
	subd.style.overflow = 'scroll' ;
	subd.style.backgroundColor = 'white' ;
	
	var f = document.createElement ( 'form' ) ;
	f.id = 'white_wizard_form' ;

	for ( var cnt = 0 ; cnt < search_results.length ; cnt++ ) {
		var s = search_results[cnt].split("\t") ;
		if ( s.length < 7 ) continue ;
		var label_text = " [" + s[5] + ":" + s[1] + "-" + s[2] + " (" + s[0] + ")]" ;
		var i = document.createElement ( 'i' ) ;
		i.appendChild ( document.createTextNode ( s[6] ) ) ;
		var label = document.createElement ( 'label' ) ;
		label.appendChild ( i ) ;
		label.appendChild ( document.createTextNode ( label_text ) ) ;
		var cb = document.createElement ( 'input' ) ;
		cb.id = "search_cb_" + cnt ;
		cb.type = 'checkbox' ;
		if ( is_MSIE() ) label.setAttribute ( 'htmlFor' , cb.id ) ;
		else label.setAttribute ( 'for' , cb.id ) ;

		f.appendChild ( cb ) ;
		f.appendChild ( label ) ;
		f.appendChild ( document.createElement ( 'br' ) ) ;
	}
	
	var b_ok = document.createElement ( 'input' ) ;
	b_ok.type = 'button' ;
	b_ok.value = i18n['ok'] ;
	b_ok.onclick = white_wizard_ok ;
	
	var b_cancel = document.createElement ( 'input' ) ;
	b_cancel.type = 'button' ;
	b_cancel.value = i18n['cancel'] ;
	b_cancel.onclick = white_wizard_cancel ;
	
	subd.appendChild ( f ) ;
	setOpacity ( subd , 100 ) ;
	setOpacity ( d , 100 ) ;
	
	d.appendChild ( subd ) ;
	d.appendChild ( b_ok ) ;
	d.appendChild ( b_cancel ) ;
	darken.appendChild ( d ) ;
}

// Event handler : OK button on search display has been clicked.
function white_wizard_ok () {
	var f = document.getElementById ( 'white_wizard_form' ) ;
	var chr = '' ;
	var from = 9999999999 ;
	var to = 0 ;
	for ( var i = 0 ; i < f.childNodes.length ; i++ ) {
		if ( f.childNodes[i].type != 'checkbox' ) continue ;
		if ( !f.childNodes[i].checked ) continue ;
		var num = parseInt ( f.childNodes[i].id.split("_").pop() ) ;
		var s = search_results[num].split ( "\t" ) ;
		if ( chr == '' ) chr = s[5] ;
		if ( s[5] != chr ) continue ; // First chromosome selected gets the pie
		if ( from > s[1] ) from = parseInt ( s[1] ) ;
		if ( to < s[2] ) to = parseInt ( s[2] ) ;
	}
	white_wizard_cancel () ;
	if ( to == 0 ) return ; // Nothing selected
	var sel = document.getElementById("chr_list");
	var chr_len = 0 ;
	for ( var i = 0 ; i < sel.options.length ; i++ ) {
		if ( sel.options[i].value != chr ) continue ;
		sel.selectedIndex = i ;
		chr_len = chromosome_length[sel.options[i].value] ;
		break ;
	}
	
	if ( chr_len > 0 ) {
		var w = Math.floor ( ( to - from + 1 ) / 6 ) ;
		from -= w ;
		to += w ;
		if ( from < 1 ) from = 1 ;
		if ( to > chr_len ) to = chr_len ;
	}
	
	document.getElementById('chr_from').value = from ;
	document.getElementById('chr_to').value = to ;
	update_image () ;
}

// Event handler : OCancel button on search display has been clicked.
function white_wizard_cancel () {
	var darken = document.getElementById('darkenScreenObject') ;
	var d = document.getElementById ( 'white_wizard' ) ;
	darken.removeChild ( d ) ;
	fade ( 'darkenScreenObject' , 50 , -1 ) ;
}

// Parses the current URL for a parameter and returns the value
function get_url_parameter( name ) {
	name = name.replace(/[\[]/,"\\\[").replace(/[\]]/,"\\\]");
	var regexS = "[\\?&]"+name+"=([^&#]*)";
	var regex = new RegExp ( regexS );
	var results = regex.exec ( window.location.href );
	if ( results == null ) return "";
	return results[1];
}

function show_second_image () {
	document.getElementById("display_second_track").checked = 1 ;
	document.getElementById('squeeze_tracks').disabled = false ;
	document.getElementById('second_image').style.display = 'block' ;
	document.getElementById('second_image_name').style.display = 'block' ;
	removeChildrenFromNode ( document.getElementById('second_image_name') ) ;
	document.getElementById('second_image_name').appendChild ( document.createTextNode ( second_track_lanes ) ) ;
}

// Initializes everything!
function init () {
	document.onmousedown = selectmouse ;
	document.onmouseup = release_drag ;
	setOpacity ( document.getElementById('lanes_container') , 90 ) ;
	if ( window.initialize_organism ) initialize_organism () ;
	initalize_chromosomes () ;
	initialize_display () ;
	show_lanes () ;
	
	document.getElementById('span_sanger_bars').style.display = sanger_layout ? 'inline' : 'none' ;
	//document.getElementById('img_container').style.width = img_width + 'px' ;

	// var sel = document.getElementById("image_width");
	// for ( var i = 0 ; i < sel.options.length ; i++ ) {
	// 	if ( sel.options[i].value >= img_width ) sel.selectedIndex = i ;
	// }
	
	if ( use_init ) {
        if ( init_max_window )
            max_window = init_max_window;

		document.getElementById('chr_from').value = init_from ;
		document.getElementById('chr_to').value = init_to ;
		
		var init_second_image = get_url_parameter ( 'second_image' ) ;
		var init_squeeze_tracks = get_url_parameter ( 'squeeze_tracks' ) ;
		
		if ( init_second_image != '' ) {
			second_track_lanes = init_second_image ;
			show_second_image () ;
			if ( init_squeeze_tracks ) document.getElementById("squeeze_tracks").checked = 1 ;
		}
		
		var cidx = 0 ;
		for ( var i = 0 ; i < chromosomes.length ; i++ ) {
			if ( chromosomes[i] != init_chr ) continue ;
			cidx = i ;
			break ;
		}
		document.getElementById("chr_list").selectedIndex = cidx ;
		
		if ( init_mode == 'pileup' ) {
			document.getElementById("display_mode_pileup").checked = true ;
			document.getElementById("display_mode_pileup").selected = true ;
			display_mode = 'pileup' ;
			document.getElementById('indel_zoom').disabled = true ;
		}
		if ( init_mode == 'paired_pileup' ) {
			document.getElementById("display_mode_paired_pileup").checked = true ;
			document.getElementById("display_mode_paired_pileup").selected = true ;
			display_mode = 'paired_pileup' ;
			document.getElementById('indel_zoom').disabled = true ;
		}
		if ( init_mode == 'coverage' ) {
			document.getElementById("display_mode_coverage").checked = true ;
			document.getElementById("display_mode_coverage").selected = true ;
			display_mode = 'coverage' ;
			document.getElementById('indel_zoom').disabled = true ;
		}

		var lanes = init_lane.split ( ',' ) ;
		if ( lanes.length == 1 ) {
			var id = 'id4label_' + init_lane ;
			document.getElementById(id).checked = true ;
			activate_lane ( init_lane ) ;
		} else {
			toggle_lanes_radio_checkbox () ;
			for ( var i = 0 ; i < lanes.length ; i++ ) {
				var id = 'id4label_' + lanes[i] ;
				document.getElementById(id).checked = true ;
			}
			update_image () ;
		}
		return ;
	}
	
	full_chromosome_rezise () ;
}

// Event handler : Double-click on feature in annotation.
function annotation_clicked ( event ) {
	if ( is_loading > 0 ) return false ;

	grayOut ( 1 , {'opacity':'0'} ) ;
	setTimeout ( "fade('darkenScreenObject',0,1);" , 1 ) ;

	var sel = document.getElementById("chr_list");
	var oldpos = document.getElementById('img_container').style.position ;
	document.getElementById('img_container').style.position = '' ;
    var pos_x = event.offsetX?(event.offsetX):event.pageX-document.getElementById("annotation_image").offsetLeft;
    var pos_y = event.offsetY?(event.offsetY):event.pageY-document.getElementById("annotation_image").offsetTop;
	document.getElementById('img_container').style.position = oldpos ;

	var url = cgi_path + '/get_data.pl' ;
	url += '?from=' + cur_from ;
	url += '&to=' + cur_to ;
	url += '&chr=' + sel.options[sel.selectedIndex].value ;
	url += '&output=url' ;
	url += '&width=' + img_width ;
	url += '&view=annotation' ;
	url += '&x=' + pos_x ;
	url += '&y=' + pos_y ;

	var request =  new XMLHttpRequest();
	request.open("GET", url, false);
	request.send(null);
	
	if ( request.responseText == '' ) {
		grayOut ( 0 ) ;
		return ;
	}
	
	update_viewport_size () ;
	var darken = document.getElementById('darkenScreenObject') ;
	var centerx = Math.floor ( viewportwidth / 2 ) ;
	var centery = Math.floor ( viewportheight / 2 ) ;
	var w2 = 200 ; // Half the width of the white box
	var h2 = 200 ; // Half the height of the white box
	var d = document.createElement ( 'div' ) ;
	d.id = 'white_wizard' ;
	d.style.zindex = 51 ;
	d.style.position = 'absolute' ;
	d.style.left = ( centerx - w2 ) + "px" ;
	d.style.top = ( centery - h2 ) + "px" ;
	d.style.width = ( w2 * 2 ) + "px" ;
	d.style.height = ( h2 * 2 ) + "px" ;
	d.style.backgroundColor = 'white' ;
	
	var lines = request.responseText.split ( "\n" ) ;
	var ul = document.createElement ( 'ul' ) ;
	ul.style.height = ( h2 * 2 - 60 ) + "px" ;
	ul.style.overflow = 'scroll' ;
	for ( var i = 0 ; i < lines.length ; i++ ) {
		var li = document.createElement ( 'li' ) ;
		li.innerHTML = lines[i] ;
		ul.appendChild ( li ) ;
	}
	
	var bcontainer = document.createElement ( 'div' ) ;
	bcontainer.style.textAlign = 'center' ;
	var b_ok = document.createElement ( 'input' ) ;
	b_ok.type = 'button' ;
	b_ok.value = i18n['ok'] ;
	b_ok.onclick = white_wizard_cancel ;

	bcontainer.appendChild ( b_ok ) ;
	d.appendChild ( ul ) ;
	d.appendChild ( bcontainer ) ;
	darken.appendChild ( d ) ;
}

//_________________________________________________________________________________________________
// Image dragging functions
//_________________________________________________________________________________________________


function movemouse(e)
{
  if (isdrag)
  {
	var newx = nn6 ? tx + e.clientX - drag_x : tx + event.clientX - drag_x ;
	var new_from = get_dragged_image_new_from ( newx ) ;
	if ( new_from < 1 ) return ;
	var cur_diff = cur_to - cur_from ;
    var max_to = get_max_to(0);
	var new_to = new_from + cur_diff ;
	if ( new_to >= max_to ) return ;
    document.getElementById('main_image').style.left = newx + "px" ;
    document.getElementById('second_image').style.left = newx + "px" ;
    document.getElementById('annotation_image').style.left = newx + "px" ;
    document.getElementById('gc_image').style.left = newx + "px" ;
    document.getElementById('coverage_image').style.left = newx + "px" ;

	var md = document.getElementById ( 'mini_data' ) ;
	removeChildrenFromNode ( md ) ;
	md.appendChild ( document.createTextNode ( new_from + "-" + new_to ) ) ;
	md.style.display = 'inline' ;

    return false;
  }
}

function selectmouse(e) 
{
	if ( is_loading > 0 ) return ;
  var fobj       = nn6 ? e.target : event.srcElement;
  var topelement = nn6 ? "HTML" : "BODY";

  while (fobj.tagName != topelement && fobj.className != "dragme")
  {
    fobj = nn6 ? fobj.parentNode : fobj.parentElement;
  }

  if (fobj.className=="dragme")
  {
	document.getElementById('vruler').style.display = 'none' ;
    isdrag = true;
    dobj = fobj;
    tx = parseInt(dobj.style.left+0);
    ty = parseInt(dobj.style.top+0);
    drag_x = nn6 ? e.clientX : event.clientX;
    drag_y = nn6 ? e.clientY : event.clientY;
	document.getElementById('img_container').style.position = 'relative' ;
    document.onmousemove=movemouse;
    return false;
  }
}

function get_dragged_image_new_from ( x ) {
	var new_left = parseInt ( x ) ;
	var cur_diff = cur_to - cur_from ;
	return Math.floor ( cur_from - cur_diff * new_left / parseInt ( img_width ) ) ;
}

function release_drag () {
	if ( !isdrag ) return false ;
	isdrag = false ;

	var new_from = get_dragged_image_new_from ( document.getElementById('main_image').style.left )  ;
	if ( new_from < 1 ) new_from = 1 ;
	if ( new_from == cur_from ) {
		document.getElementById('img_container').style.position = '' ;
		return ;
	}
	
	var cur_diff = cur_to - cur_from ;
	var new_to = Math.floor ( new_from + cur_diff ) ;
	
	document.getElementById('chr_from').value = new_from ;
	document.getElementById('chr_to').value = new_to ;
	update_image () ;
}

function togglesidebar() {
	var state = document.getElementById('lanes').style.display ;
	if ( state == 'none' ) {
		document.getElementById('lanes').style.display = 'block' ;
	} else {
		document.getElementById('lanes').style.display = 'none' ;
	}
	return false ;
}
