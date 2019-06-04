
$(window).load(function() {
    // on dom ready
if(document.getElementById("cy")){  
    var disease_on = $("#disease_on");
    var gene_on = $("#gene_on");
    var disease_names_on = $("#show_disease_names");
    var gene_names_on = $("#show_gene_names");
    var save_photo = $("#save_photo");
    disease_on.bootstrapSwitch();
    gene_on.bootstrapSwitch();
    disease_names_on.bootstrapSwitch();
    gene_names_on.bootstrapSwitch();
    
        $("#cy").cytoscape({
            style:cytoscape.stylesheet().selector("node").css({
                content:"data(id)",
                "text-valign":"center",
                "font-size":"mapData(weight, 0, 1, 18px, 22px)",
                color:"white",
                "text-outline-width":"mapData(weight, 0, 0.5, 4px, 3px)"
            }).selector("edge").css({
                "target-arrow-shape":"none",
                "line-color":"#444401",
                width:"0.5px"
            }).selector("edge.BIOSYSTEM").css({
                "line-color":"#baba11"
            }).selector("edge.HPRD").css({
                "line-color":"mapData(weight,0.2,0.8, #00ef33, #00ef33)"
            }).selector("edge.GENE_FAMILY").css({
                "line-color":"#6600CC"
            }).selector("edge.HTRI").css({
                "line-color":"#020202",
                "target-arrow-shape":"triangle",
                "target-arrow-color":"#1f1f1f"
            }).selector("edge.disease_term").css({
                "line-color":"#0099ff"
            }).selector("edge.disease_gene").css({
                "line-color":"#0099ff"
            }).selector(":selected").css({
                "background-color":"red",
                "line-color":"red",
                "target-arrow-color":"red",
                "source-arrow-color":"red"
            }).selector("node.Reported").css({
                "background-color":"mapData(color_weight,0,1,#78cea3,#00cc66)",
                "text-outline-color":"mapData(color_weight,0,1,#78cea3,#00cc66)",
                "border-width":"0px",
                width:"mapData(weight, 0,1, 16, 85)",
                height:"mapData(weight, 0,1, 16, 85)"
            }).selector("node.Predicted").css({
                "background-color":"mapData(color_weight,0,1, #f4bb83, #FF9933)",
                "text-outline-color":"mapData(color_weight, 0, 1, #f4bb83, #FF9933)",
                "border-width":"0px",
                width:"mapData(weight, 0,1, 16, 85)",
                height:"mapData(weight, 0,1, 16, 85)"
            }).selector("node.disease").css({
                "text-outline-width":"2px",
                "text-outline-color":"mapData(weight,0,1,#a3daff,#0099ff)",
                shape:"roundrectangle",
                "background-color":"mapData(weight,0,1,#a3daff,#0099ff)",
                width:"mapData(weight, 0, 1, 18, 90)",
                height:"20",
                "font-size":"10px"
            }).selector("node.term").css({
                "text-outline-width":"3px",
                "text-outline-color":"#FF33CC",
                "background-color":"#fe2200",
                "font-size":"35px",
                shape:"circle",
                width:"25px",
                height:"25px"
            }).selector(".highlighted").css({
                "background-color":"#ff0303",
                "line-color":"#ff0303",
                "text-outline-color":"#ff0303",
                "target-arrow-color":"#ff0303",
                "source-arrow-color":"#ff0303"
            }).selector(".faded").css({
                opacity:"0"
            }),
            userZoomingEnabled:true,
            panningEnabled:true,
            userPanningEnabled:true,
            layout:{
                name:"random"
            },
            zoom:1,
            minZoom:1e-50,
            maxZoom:1e50,
            elements:[],
            ready:function() {
                window.cy = this;
                cy.elements().unselectify();
                cy.on("tap", "node", function(e) {
                    var node = e.cyTarget;
                    var neighborhood = node.neighborhood().add(node);
                    //  cy.elements().removeClass('highlighted');
                    // neighborhood.addClass('highlighted');
                    cy.elements().addClass("faded");
                    neighborhood.removeClass("faded");
                    if (node.hasClass("disease")) {
                        if (disease_names_on.prop("checked") == false) node.css("content", "data(id)");
                    }
                });
                cy.on("tap", function(e) {
                    if (e.cyTarget === cy) {
                        //cy.elements().removeClass('highlighted');
                        cy.elements().removeClass("faded");
                        if (disease_names_on.prop("checked") == false) cy.elements("node.disease").css("content", "");
                    }
                });
            }
        });
        var cy = $("#cy").cytoscape("get");

		cy.load(json1);
        cy.boxSelectionEnabled(false);
        var defaults = {
            zoomFactor:.05,
            // zoom factor per zoom tick
            zoomDelay:45,
            // how many ms between zoom ticks
            minZoom:.1,
            // min zoom level
            maxZoom:10,
            // max zoom level
            fitPadding:50,
            // padding when fitting
            panSpeed:10,
            // how many ms in between pan ticks
            panDistance:10,
            // max pan distance per tick
            panDragAreaSize:75,
            // the length of the pan drag box in which the vector for panning is calculated (bigger = finer control of pan speed and direction)
            panMinPercentSpeed:.25,
            // the slowest speed we can pan by (as a percent of panSpeed)
            panInactiveArea:8,
            // radius of inactive area in pan drag box
            panIndicatorMinOpacity:.5,
            // min opacity of pan indicator (the draggable nib); scales from this to 1.0
            autodisableForMobile:true,
            // disable the panzoom completely for mobile (since we don't really need it with gestures like pinch to zoom)
            // icon class names
            sliderHandleIcon:"fa fa-minus",
            zoomInIcon:"fa fa-plus",
            zoomOutIcon:"fa fa-minus",
            resetIcon:"fa fa-expand"
        };
        // $("#cy").cytoscapePanzoom(defaults);
        //set the initial states of checkboxes  
        disease_on.prop("checked", true);
        gene_on.prop("checked", true);
        gene_names_on.prop("checked", true);
        disease_names_on.prop("checked", false);
        var disease_nodes;
        var gene_nodes;
        var disease_edges = cy.elements("edge.disease_term");
        var disease_gene_edges = cy.elements("edge.disease_gene");
        var gene_edges = cy.elements("edge.gene");
        var term_nodes;
        disease_on.on("switchChange.bootstrapSwitch", function() {
            if (!disease_on.is(":checked")) {
                disease_nodes = cy.elements("node.disease");
                term_nodes = cy.elements("node.term");
                cy.remove(disease_nodes);
                cy.remove(term_nodes);
                layout_change(1e3);
            } else {
                cy.add(disease_nodes);
                cy.add(term_nodes);
                cy.add(disease_edges);
                cy.add(disease_gene_edges);
                layout_change(1e3);
            }
        });
        gene_on.on("switchChange.bootstrapSwitch", function() {
            if (!gene_on.is(":checked")) {
                gene_nodes = cy.elements("node.gene");
                cy.remove(gene_nodes);
                layout_change(1e3);
            } else {
                cy.add(gene_nodes);
                cy.add(gene_edges);
                cy.add(disease_gene_edges);
                layout_change(1e3);
            }
        });
        disease_names_on.on("switchChange.bootstrapSwitch", function() {
            if (disease_names_on.is(":checked")) {
                cy.elements("node.disease").css("content", "data(id)");
            } else {
                cy.elements("node.disease").css("content", "");
            }
        });
        gene_names_on.on("switchChange.bootstrapSwitch", function() {
            if (gene_names_on.is(":checked")) {
                cy.elements("node.gene").css("content", "data(id)");
            } else {
                cy.elements("node.gene").css("content", "	");
            }
        });
        var show_edges = $("#show_edges");
        show_edges.change(function() {
            if (show_edges.find("option:selected").val() == "all") {
                cy.elements("edge").show();
            }
            if (show_edges.find("option:selected").val() == "HPRD") {
                cy.elements("edge").hide();
                cy.elements("edge.HPRD").show();
            }
            if (show_edges.find("option:selected").val() == "GeneFamily") {
                cy.elements("edge").hide();
                cy.elements("edge.GENE_FAMILY").show();
            }
            if (show_edges.find("option:selected").val() == "Biosystem") {
                cy.elements("edge").hide();
                cy.elements("edge.BIOSYSTEM").show();
            }
            if (show_edges.find("option:selected").val() == "TranscriptionInteraction") {
                cy.elements("edge").hide();
                cy.elements("edge.HTRI").show();
            }
        });
		$("#undo").click(function(){
                        if (adjust_layout.find("option:selected").val() != "circle"||"grid") {
                cy.layout({
                    name:"random"
                });
            };
if (adjust_layout.find("option:selected").val() == "circle") {
                cy.layout({
                    name:"circle",
                    avoidOverlap:true
                });
            };
            if (adjust_layout.find("option:selected").val() == "grid") {
                cy.layout({
                    name:"grid",
                    avoidOverlap:false
                });
            }
		});
        $("#save_photo").click(function() {
            window.open(cy.png({
                bg:"#ffffff"
            }));
        });
        var adjust_layout = $("#adjust_layout");
        function layout_change(time) {
            if (adjust_layout.find("option:selected").val() == "force") {
                cy.layout({
                    name:"random"
                });
            };
            if (adjust_layout.find("option:selected").val() == "circle") {
                cy.layout({
                    name:"circle",
					avoidOverlap:false
                });
            };
            if (adjust_layout.find("option:selected").val() == "grid") {
                cy.layout({
                    name:"grid",
                    avoidOverlap:false
                });
            }

        }
        adjust_layout.change({
            time:1e3
        }, function(event) {
            layout_change(event.data.time);
        });
        /* setTimeout(function() {
            adjust_layout.selectpicker("val", "force");
            adjust_layout.selectpicker("refresh");
            adjust_layout.change();
        }, 1500); */
    
}
});