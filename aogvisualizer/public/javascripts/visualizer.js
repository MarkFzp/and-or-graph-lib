/*
 * Description: Visualizer of And-Or Graph with D3.js
 * Author: Feng Shi
 * Date of Creation: July 12th, 2017
 * 
 */

/////////////// Data parsing ///////////////////
// var ast, graph;
var str;
var colors = { "And-node":      "#111",
               "Or-node":       "#ccc",
               "stroke-mover":  "#F4D03F",
               "stroke-normal": "#566573" };
var w = 1200, h = 600;
var nodeRadius = 6;

var svg; // = d3.select("#chart").append("svg:svg").attr("width", w).attr("height", h);

var inFile = document.getElementById("inFile");

inFile.addEventListener('change', function () {
	if (!d3.select("#chart").select("svg").empty())
		d3.select("#chart").select("svg").remove();
	
	//svg = d3.select("#chart").append("svg:svg").attr("width", w).attr("height", h);

    svg = d3.select("#chart").append("svg:svg")
            .attr("width", w)
            .attr("height", h)
            .append("g")
            .call(d3.behavior.zoom().scaleExtent([0.1, 8]).on("zoom", zoom))
            .append("g");

    svg.append("rect")
        .attr("class", "overlay")
        .attr("width", w)
        .attr("height", h);

    function zoom() {
        svg.attr("transform", "translate(" + d3.event.translate + ")scale(" + d3.event.scale + ")");
    }

	if (!inFile.files.length) {
		alert("Please select a file!");
		return;
	}

	var fi = inFile.files[0];
	var reader = new FileReader();
	str = "";
	
	var textType = /text.*/;
	if (fi.type.match(textType)) {
		reader.onload = function() {
			str = reader.result;
			createAoG(str);
		};
		reader.readAsText(fi);
	}
});

function createAoG (str) {
	// parse the dot file
	// var ast = DotParser.parse(str);
	// var graph = new DotGraph(ast);
	// graph.walk();
    //
	// // translate dot file into Json compatible format
	// var gr = new JGraph(graph.nodes, graph.edges);
	var gr = dot2json(str);

    var weights = getWeights(gr);

	gr.root.fixed = true;
	gr.root.x = w / 2;
	gr.root.y = h / 10;

    // build the arrow marker.
    svg.append("svg:defs").selectAll("marker")
        .data(["arrow0", "arrow1"])      // Different link/path types can be defined here
        .enter().append("svg:marker")    // This section adds in the arrows
        .attr("id", String)
        .attr("viewBox", "0 -5 10 10")
        .attr("refX", 20)
        .attr("refY", -0.5)
        .attr("markerWidth", 3)
        .attr("markerHeight", 3)
        .attr("orient", "auto")
        .append("svg:path")
        .attr("d", "M0,-5L10,0L0,5");


	var nodes = flattenD3(gr.root),
	// var nodes = d3.layout.tree(gr);
		links = d3.layout.tree().links(gr.nodes);
	//////
	var force = d3.layout.force()
                .gravity(.2)
                .charge(-1500)
                .size([ w, h ])
				.linkDistance(30);

	var link = svg.selectAll("line.link").data(links).enter()
				.append("svg:line")
				.attr("class", "link")
                .style("display", function(d) {
                    if (d.source.AO === 'X') {
                        return "none";
                    } // fake root;
                });

	var node = svg.selectAll("g.node").data(nodes).enter()
			    .append("svg:g")
			    .attr("class", "node")
			    .attr("type", function(d) { return d.AO;   })
			    .attr("name", function(d) { return d.name; })
			
	node.append("svg:circle")
		.attr("r", nodeRadius)
		.style("fill", function(d) {
			if (d.AO === 'A') {
                return colors["And-node"];
            } // d.color;
		})
        .style("display", function(d) {
            if (d.AO === 'X') {
                return "none";
            } // fake root;
        })
		.call(force.drag);

	node.append("svg:text")
        .style("display", function(d) {
            if (d.AO === 'X') {
                return "none";
            } // fake root;
        })
		.style("pointer-events", "none")
		.text(function (d) { return d.label; })
		.attr("text-anchor", "middle")
		.attr("x", 10)
		.attr("y", -18)
		.style("font-size", 15);

	/// add an arrow to each link
	d3.selectAll("line.link").attr("marker-end", "url(#arrow0)");

    force.nodes(nodes)
		.links(links)
		.on("tick", tick)
		.start();

    function tick(e) {
        var ky = 1.2 * e.alpha;
        links.forEach(function(d) {
            d.target.y += (d.target.depth * 120 - d.target.y) * ky;
        });

        var nds = nodes.filter(function (t) { return t.AO === "A"; });
        nds.forEach(function (t) {
            // each and-node may have children with spatial position
            if (t.hasOwnProperty("children") && (t.children.length > 1))
                spatialSort(t.children);
        });

        d3.selectAll("g.node")
            .attr("transform", function (d) {
                return "translate(" + d.x + "," + d.y + ")";
            });

        d3.selectAll("line.link")
            .attr("x1", function(d) { return d.source.x; })
            .attr("y1", function(d) { return d.source.y; })
            .attr("x2", function(d) { return d.target.x; })
            .attr("y2", function(d) { return d.target.y; });
    }

	function spatialSort(nodes) {
        var xs = [];
        // 1. sort the order of x position
        for (var key in nodes) { xs.push(nodes[key].x); }
        xs.sort(function (a, b) { return a - b; });
        // 2. sort the nodes with their defined spatial position
        nodes.sort(function (a, b) { return a.pos - b.pos; });
        // 3. assign the corresponding position to each node
        for (var key in nodes) { nodes[key].x = xs[key]; }
    }

	var lks = [], nds = [];

    // node.on("click", function(d) {
	// 	/// Avoid dragged event
	// 	if (d3.event.defaultPrevented) return;
	//
	// 	/// clean previously selected nodes, links and arrows
	// 	if (nds.length !== 0) {
	// 		nds.selectAll("circle").style("fill", function (n) {
     //            /** @namespace n.AO */
     //            if (n.AO === "A") {
     //                return colors["And-node"];
     //            } // n.color;
	// 		});
	//
	// 		lks.style("stroke", function (l) {
	// 			return l.stroke;
	// 		}).attr("marker-end", "url(#arrow0)");
    //
	// 	}
	//
	// 	var pg  = sample(d);
	//
	// 	nds = node.filter(function (n) {
	// 					return (pg.indexOf(n) !== -1);
	// 				});
	//
	// 	nds.selectAll("circle").style("fill", "yellow");
	//
	// 	lks = link.filter(function (l) {
	// 		return ((pg.indexOf(l.source) !== -1) && (pg.indexOf(l.target) !== -1));
	// 	});
	//
	// 	lks.style("stroke", "red").attr("marker-end", "url(#arrow1)");
	//
	// });

    var tooltip = d3.select("body").append("div")
        .attr("class", "tooltip")
        .style("opacity", 0);

	link.on("mouseover", function(d) {
		d3.select(this).attr("marker-end", "url(#arrow1)")
            .style("stroke", "red");

		tooltip.transition()
			.duration(200)
			.style("display", "")
            .style("opacity", .9);
        tooltip.html("weight: " + weights[d.source.label+"to"+d.target.label])
            .style("left", (d3.event.pageX + 12) + "px")
            .style("top", (d3.event.pageY - 28) + "px");
    })
    .on("mouseout", function(d) {
		d3.select(this).attr("marker-end", "url(#arrow0)")
			.style("stroke", "#66ccff");

        tooltip.transition()
			.duration(200)
			.style("display", "none");
    });

	node.on("mouseover", function (d) {
            d3.select(this).select("circle").transition()
                .duration(500)
                .attr("r", 3*nodeRadius)
                .style("stroke-width", "4px")
                .style("stroke", colors["stroke-mover"]);
            d3.select(this).select("text").transition()
                .duration(500)
                .style("stroke-width", ".5px")
                .style("font", "26.8px serif")
                .style("opacity", 1);

            tooltip.transition()
        	    .duration(200)
        	    .style("display", "")
        	    .style("opacity", .9);
        	tooltip.html(function() {
        			if(d.AO === "A")
        				return "And-node: " + d.label;
        			return "Or-node: " + d.label;
        		})
        		.style("background", "AAAAAA")
            	.style("left", (d3.event.pageX + 12) + "px")
            	.style("top", (d3.event.pageY - 28) + "px");
        })
        .on("mouseout", function () {
            d3.select(this).select("circle").transition()
                .duration(500)
                .attr("r", nodeRadius)
                .style("stroke-width", "1px")
                .style("stroke", colors["stroke-normal"]);
            d3.select(this).select("text").transition()
                .duration(500)
                .style("font", "16px Times New Roma")
                .style("opacity", 0.8);
            tooltip.transition()
            	.duration(200)
				.style("display", "none");
        });
}

function flattenD3(root) {
	var nodes = [];
	function recurse(node, depth, seen) {
		if (seen.indexOf(node) > -1) {
			return;
		}
		if (node.children) {
			node.children.forEach(function(child) {
				recurse(child, depth + 1, seen.concat([node]));
			});
		}
		// if (node.depth < depth) node.depth = depth; // as latest as possible
		node.depth = depth; // as earliest as possible
		nodes.push(node);
	}
	recurse(root, 1, []);
	return nodes;
}

function sample(start) {
	var pg = [];
	function recurse (node) {
		if(node.children) {
			var idx;
			switch (node.AO) {
			case "A": 
				for (idx in node.children) { recurse(node.children[idx]); }
				break;
			case "O": 
				idx = Math.ceil(Math.random() * node.children.length) - 1;
				recurse(node.children[idx]);
				break;
			default:
				break;
			}
		}
		pg.push(node);
	}
	recurse(start);
	return pg;
}

function getWeights(graph){
    var weight = {};
    for(var i = 0; i<graph.edges.length; i++) {
        var src = graph.edges[i].source.label;
        var tgt = graph.edges[i].target.label;
        var wgt = graph.edges[i].weight;

        if(weight[src+"to"+tgt] === undefined)
            weight[src+"to"+tgt] = [wgt];
        else
            weight[src+"to"+tgt].push(wgt);
    }

    return weight;
}

////////Create a legend //////////
var legdata = ["And-node", "Or-node"];
var legvis  = d3.select("#footer").append("svg:svg").attr("width", 150).attr("height", 100);
var legend = legvis.append("g")
			.attr("class", "legend")
			.attr("x", 20) // w - 65)
			.attr("y", h - 100)
			.attr("height", 100)
			.attr("width", 100);

legend.selectAll('g').data(legdata)
	.enter()
	.append('g')
	.each(function (d, i) {
		var g = d3.select(this);
			g.append("circle")
			 .attr("r", nodeRadius)
			 .attr("cx", 15) // w - 65)
			 .attr("cy", 15 + i * 25)
			 .style("fill", colors[legdata[i]])
			 .style("stroke", "black")
			 .style("stroke-width", 1);
			g.append("text")
			 .attr("x", 40) // w - 50)
			 .attr("y", 20 + i * 25)
			 .attr("height", 30)
			 .attr("width", 50)
			 .style("fill", "black")
			 .style("font-family", "Arial")
			 .style("font-size", 12)
			 .text(legdata[i]);
	});
///////////////////////////////////

