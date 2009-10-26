var xpos, ypos;

jQuery.noConflict();

jQuery(document).ready(init);

function init()
{
    jQuery().mousemove(mouse_moved);
    // jQuery("#gene").keyup(function(e){e.preventDefault(); ajax_search(); }); 
}

function ajax_search()
{ 
    jQuery.post("snps.pl?ajax=1&gene=" + jQuery("#gene").val() , function(data){
            var html = '';
            if (data.length>0)
            {
                var list = eval(data);
                for (var i=0; i<list.length; i++)
                    html += "<span style='cursor:pointer' onclick=\"set_value('#gene','" + list[i] +
                        "');hide_element('#gene_results');\">" +list[i]+ "</span><br />";
            }
            jQuery("#gene_results").html(html); 
            jQuery("#gene_results").show(); 
            });
}

function set_value(element,value)
{
    jQuery(element).attr('value',value);
}

function mouse_moved(event)
{
    if ( event.pageX || event.pageY )
    {
        xpos = event.pageX;
        ypos = event.pageY;
    }
}

function fadeout_element(name,ms)
{
    if ( ms==null ) ms = 1000;
    setTimeout(function(){ jQuery(name).fadeOut('fast') }, ms);
}
function show_element(name)
{
    jQuery(name).show();
}
function hide_element(name)
{
    jQuery(name).hide();
}
function toggle_element(name)
{
    jQuery(name).toggle();
}
function toggle_element_here(name)
{
    var height = jQuery(name).height();
    var width  = jQuery(name).width();
    jQuery(name).css("left", xpos-width/2 + "px");
    jQuery(name).css("top", ypos-height + "px");
    jQuery(name).toggle();
}

var last_toggled=0;
function toggle_one_element(name)
{
    jQuery(name).toggle();
    if ( last_toggled && last_toggled!=name )
        jQuery(last_toggled).hide();
    last_toggled = name;
}

function show_mouse_info(strain)
{
    jQuery('#mouse_info').html(strain);
}

function check_all(prefix,value)
{
    jQuery("input[name^='"+prefix+"']").attr('checked',value);
}

