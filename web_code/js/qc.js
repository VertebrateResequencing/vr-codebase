jQuery(document).ready(function() {

        //jQuery("a").click(function() { alert("jQuery OK"); }); // Debug

        jQuery("#resetRadios").click(function() {
            // reset all the radio buttons

            jQuery("input:radio").attr("checked", false);
            });

        jQuery("input:radio[id*='toggle']").click(function(){
            // toggle all radios of this lane status
    
            status = jQuery(this).attr('id');
            status = status.replace("toggle_","");

            state = jQuery(this).is(':checked');
            jQuery(this).attr('checked',false);

            // this failes 2nd time round, does not reset checkbox
            jQuery("input:radio").each(function(){
                if (jQuery(this).attr('value') == status ) {
                    jQuery(this).attr('checked',state);
                }
            });
        });

        jQuery("table[id*='srt']").tablesorter();

});

