var str = 'digraph g {'
		+ '10[label="10_10", name="10_10", shape="oval", color="#8888ff", STC="S", style="bold", AO="O"];'
		+ '3[label="3_3", name="3_3", shape="oval", color="#8888ff", STC="S", style="bold", AO="A"];'
		+ '3 -> 7;'
		+ '2[label="2_2", name="2_2", shape="oval", color="#8888ff", STC="S", style="filled", AO="A"];'
		+ '2 -> 3;'
		+ '2 -> 4;'
		+ '1[label="1_1", name="1_1", shape="oval", color="#8888ff", STC="S", style="bold", AO="O"];'
		+ '1 -> 8;'
		+ '0[label="0_0", name="0_0", shape="oval", color="#8888ff", STC="S", style="filled", AO="A"];'
		+ '0 -> 1;'
		+ '0 -> 2;'
		+ '4[label="4_4", name="4_4", shape="oval", color="#8888ff", STC="S", style="filled", AO="O"];'
		+ '4 -> 5;'
		+ '4 -> 6;'
		+ '5[label="5_5", name="5_5", shape="oval", color="#8888ff", STC="S", style="bold", AO="O"];'
		+ '6[label="6_6", name="6_6", shape="oval", color="#8888ff", STC="S", style="filled", AO="O"];'
		+ '7[label="7_7", name="7_7", shape="oval", color="#8888ff", STC="S", style="bold", AO="O"];'
		+ '8[label="8_8", name="8_8", shape="oval", color="#8888ff", STC="S", style="filled", AO="O"];'
		+ '1 -> 9;'
		+ '1 -> 10;'
		+ '9[label="9_9", name="9_9", shape="oval", color="#8888ff", STC="S", style="bold", AO="O"];'
		+ '}';
ast = DotParser.parse(str);
graph = new DotGraph(ast);
graph.walk(); 

var gr = new jGraph(graph.nodes, graph.edges);

console.log(gr.nodes);
console.log(gr.edges);

console.log(gr.root);

var w = 680, h = 800;
var nodeRadius = 10;

//////
force = d3.layout.force().gravity(.2).charge(-1500).size([ w, h ])
		.linkDistance(30);

var svg = d3.select("#chart").append("svg:svg").attr("width", w).attr("height",
		h);

svg.append("svg:rect").attr("width", w).attr("height", h);

var marker = d3.select("svg").append('defs')
				.append('marker')
				.attr("id", "Triangle")
				.attr("refX", 12)
				.attr("refY", 6)
				.attr("markerUnits", 'userSpaceOnUse')
				.attr("markerWidth", 12)
				.attr("markerHeight", 18)
				.attr("orient", 'auto')
				.append('path')
				.attr("d", 'M 0 0 12 6 0 12 3 6');

// d3.json("flare.json", function(root) {
var nodes = flattenD3(gr.root), links = d3.layout.tree().links(gr.nodes);

gr.root.fixed = true;
gr.root.x = w / 2;
gr.root.y = 150;

var link = svg.selectAll("line").data(links).enter().insert("svg:line");

var node = svg.selectAll("circle.node").data(nodes).enter()
		.append("svg:circle")
		.attr("type", function(d) { return d.AO;   })
		.attr("name", function(d) { return d.name; })
		.attr("class", "node")
		.attr("r", nodeRadius)
		.style("fill", function(d) {
			if (d.AO === 'A')
				return d.color;
		})
		.call(force.drag);

node.append("svg:text")
	.text(function (d) { return d.label; })
	.attr("text-anchor", "")
	// .style("pointer-events", "none")
	.style("fill", "#555")
	.style("font-family", "Arial")
	.style("font-size", 12);

d3.selectAll("line").attr("marker-end", "url(#Triangle)");

force.nodes(nodes).links(links).start();

force.on("tick", function(e) {
	var ky = 1.2 * e.alpha;
	links.forEach(function(d, i) {
		d.target.y += (d.target.depth * 120 - d.target.y) * ky;
	});

	node.attr("cx", function(d) { return d.x; })
		.attr("cy", function(d) { return d.y; });

	link.attr("x1", function(d) { return d.source.x; })
		.attr("y1", function(d) { return d.source.y; })
		.attr("x2", function(d) { return d.target.x; })
		.attr("y2", function(d) { return d.target.y; });
	
});

var lks = [], nds = [];
node.on("click", function(d) {
//	if (d3.event.defaultPrevented) return; // dragged
//	node.filter(function(n) {
//		return n !== d;
//	}).classed("selected", false).attr("r", nodeRadius).style("fill",
//			function(d) {
//				if (d.AO === 'A')
//					return d.color;
//			});
	
	if (nds.length !== 0) {
		nds.style("fill", function (n) {
		if (n.AO == 'A') return n.color;
		});
		
		lks.style("stroke", function (l) {
			return l.stroke;
		});
	}
	
	// link.classed("selected", false);

	var pg  = sample(d);
	
	nds = node.filter(function (n) {
					return (pg.indexOf(n) != -1);
				});
	
	nds.style("fill", "red");
	
	lks = link.filter(function (l) {
		return ((pg.indexOf(l.source) != -1) && (pg.indexOf(l.target) != -1));
	})
	
	lks.style("stroke", "red");
	
	console.log(nds);
});

function flattenD3(root) {
	var nodes = [];
	function recurse(node, depth) {
		if (node.children) {
			node.children.forEach(function(child) {
				recurse(child, depth + 1);
			});
		}
		// if (node.depth < depth) node.depth = depth; // as latest as possible
		node.depth = depth; // as earliest as possible
		nodes.push(node);
	}
	recurse(root, 1);
	return nodes;
}

function sample(start) {
	var pg = [];
	function recurse (node) {
		if(node.children) {
			switch (node.AO) {
			case "A": 
				for (var idx in node.children) { recurse(node.children[idx]); }
				break;
			case "O": 
				var idx = Math.ceil(Math.random() * (node.children.length * 1.0)) - 1;
				recurse(node.children[idx]);
				break;
			default:
				console.log("Sample Algo: unknown option node.AO: " + node.AO);
				break;
			}
		}
		pg.push(node);
	}
	recurse(start);
	return pg;
}
