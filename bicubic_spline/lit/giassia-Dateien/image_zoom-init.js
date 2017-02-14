jQuery(document).ready(function( $ ){
    var options = IZ.options;
    $(".zoooom, .zoooom img").image_zoom(options);

    // WooCommerce category pages
    if ( IZ.woo_categories == '1' ) {
        $(".tax-product_cat .products img").image_zoom(options);
    }

    // Show zoom for lazy_load images
    if ( IZ.lazy_load == '1' ) {
        if (typeof $.unveil === "function") { 
            $("img.unveil").unveil(0, function() {
                $(this).load(function() {
                    $("img.zoooom, .zoooom img").image_zoom(options);
                });
            });
        }
    }


    // Resize the zoom windows when resizing the page
    $(window).bind('resize', function(e) {
        window.resizeEvt;
        $(window).resize(function() {
            clearTimeout(window.resizeEvt);
            window.resizeEvt = setTimeout(function() {
                $(".zoomContainer").remove();
                $(".zoooom, .zoooom img, .attachment-shop_single").image_zoom(options);
            }, 500);
        });
    });

    
    // Show zoom on the WooCommerce gallery
    if ( IZ.with_woocommerce == '1' ) { 
    $(".attachment-shop_single").image_zoom(options);

    $("a[data-rel^='zoomImage']").each(function(index){
        $(this).click(function(event){
            // If there are more than one WooCommerce gallery, exchange the thumbnail with the closest .attachment-shop_single
            var obj1 = $(".attachment-shop_single");
            if ( obj1.length > 1 ) {
                var obj1 = $(this).closest('.images').find( $(".attachment-shop_single") );
            }
            var obj2 = $(this).find("img");

            event.preventDefault();

            if ( obj2.hasClass('attachment-shop_single') === false ) {

                // Remove the srcset and sizes
                obj1.removeAttr('srcset').removeAttr('sizes');
                obj2.removeAttr('srcset').removeAttr('sizes');

                var thumb_src = obj2.attr('src');

                // Exchange the attributes
                $.each(['alt', 'title'], function(key,attr) {
                    var temp;
                    if ( obj1.attr( attr ) ) temp = obj1.attr( attr ); 
                    if ( obj2.attr( attr ) ) {
                        obj1.attr(attr, obj2.attr(attr) );
                    } else {
                        obj1.removeAttr( attr );
                    }
                    if ( temp && temp.length > 0 ) {
                        obj2.attr(attr, temp);
                    } else {
                        obj2.removeAttr( attr );
                    }
                });

                // Exchange the link sources
                var temp;
                temp = obj2.parent().attr('href');
                obj2.parent().attr('href', obj1.parent().attr('href'));
                obj1.parent().attr('href', temp );

                // Set the obj1.src = the link source
                obj1.attr('src', temp ); 

                // Set the obj2.src = data-thumbnail-src
                if ( obj1.data('thumbnail-src') ) {
                    obj2.attr( 'src', obj1.attr('data-thumbnail-src'));
                }

                // Set the obj1.data-thumbnail-src
                obj1.attr('data-thumbnail-src', thumb_src ); 

                // Remove the old zoom and reactive the new zoom
                $(".zoomContainer").remove();
                $(".attachment-shop_single").image_zoom(options);
            }

            });
        });
    }

});

