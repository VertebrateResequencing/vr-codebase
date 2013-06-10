// QC Grind jQuery

jQuery(document).ready(function() {

    // sorter for table with class sortable
    jQuery("table.sortable").tablesorter();

    // study views - reset all radio buttons and reload form
    jQuery("#resetRadios").click(function() {

        jQuery("input:radio").attr("checked", false);
        jQuery("form[name='qc_status']").submit();
        });

    // study views - toggle each lane status radio set
    jQuery("input:radio.togglestatus").click(function(){

        status = jQuery(this).attr('id');           // eg passed
        checked = jQuery(this).is(':checked');        // true/false

        jQuery("input:radio[value='" + status + "']").each(function(){
            // fails 2nd time round, browser does not reset checked state
            jQuery(this).attr('checked',checked);
        });

    });

    // AJAX Project search
    jQuery("#projFind").keyup(function() {

        var query=jQuery("#projFind").val();

        if (query == '') {
            jQuery("#projHint").html('');
            return;
        }

        jQuery.ajax({
                type: 'GET',
                url: '/cgi-bin/teams/team145/qc_grind/project_hint.pl',
                data: {"q": query},
                success: function(data) {
                    jQuery("#projHint").html(data);
                }
            });
    });

});

