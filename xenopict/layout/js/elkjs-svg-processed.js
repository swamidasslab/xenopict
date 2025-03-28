"use strict";

// XML Classes


class Node {
  constructor() {
    this.INDENTATION = 2;
  }

  indent(indent) {
    return " ".repeat(this.INDENTATION * indent);
  }
}

class Text extends Node {
  constructor(content) {
    super();
    this.content = this._stripCommonIndentation(content);
  }

  _stripCommonIndentation(content) {
    var count = 0;
    var firstChar = true;
    for (const char of content) {
      if (firstChar && char == "\n") {continue}
      if (char == " ") {count += 1}
      else if (char == "\t") {count += this.INDENTATION}
      else {break}
      firstChar = false;
    }
    return content.split("\n").filter(line => line.trim()).map(line => line.slice(count)).join("\n");
  }

  render(indent=0) {
    return this.content.split("\n").map(line => this.indent(indent) + line).join("\n");
  }
}

class Cdata extends Node {
  constructor(child) {
    super();
    this.child = child;
  }

  render(indent=0) {
    var out = [];
    out.push(`${this.indent(indent)}<![CDATA[\n`);
    out.push(this.child.render(indent + 1));
    out.push(`\n${this.indent(indent)}]]>`);
    return out.join("");
  }
}

class Xml extends Node {
  constructor(tagName, attributes={}, children=[]) {
    super();
    this.tagName = tagName;
    this.attributes = attributes;
    this.children = children;
  }

  render(indent=0) {
    var attributesString = Object.keys(this.attributes).map(key => `${key}="${this.attributes[key]}"`).join(" ")
    var numChildren = this.children.length;
    var out = [];

    out.push(`${this.indent(indent)}<${this.tagName}${attributesString? " " + attributesString: ""}`);
    if (this.children.length != 0) {
      out.push(">");
      out.push("\n");
      for (const child of this.children) {
        out.push(`${child.render(indent + 1)}\n`);
      }
      out.push(this.indent(indent));
      out.push(`</${this.tagName}>`);
    }
    else {
      out.push(` />`);
    }

    return out.join("");
  }
}

if (false) {
  console.log(
    new Xml("svg", {"version": "1.1", "xmlns": "http://www.w3.org/2000/svg", "width": 104, "height": 104}, [
      new Xml("defs", {}, [
        new Xml("style", {"type": "text/css"}, [
          new Cdata(
            new Text(`
              rect {
                opacity: 0.8;
                fill: #6094CC;
                stroke-width: 1;
                stroke: #222222;
              }

              rect.port {
                opacity: 1;
                fill: #326CB2;
              }
            `)
          )
        ])
      ]),
      new Xml("g", {"transform": "translate(0,0)"}, [
        new Xml("polyline", {"points": "42,27 62,27", "id": "e1", "class": "edge"}),
        new Xml("polyline", {"points": "42,37 52,37 52,77 62,77", "id": "e2", "class": "edge"}),
        new Xml("rect", {"id": "n1", "class": "node", "x": "12", "y": "17", "width": "30", "height": "30"}),
        new Xml("rect", {"id": "n2", "class": "node", "x": "62", "y": "12", "width": "30", "height": "30"}),
        new Xml("rect", {"id": "n3", "class": "node", "x": "62", "y": "62", "width": "30", "height": "30"}),
      ])
    ]).render()
  )
}

window.elkSvgHelpers = {
  Xml,
  Text,
  Cdata
};


// ELK SVG


// XML helpers are already in global scope

function Renderer() {
  // configuration
  this._style = `
    rect {
      opacity: 0.8;
      fill: #6094CC;
      stroke-width: 1;
      stroke: #222222;
    }
    rect.port {
      opacity: 1;
      fill: #326CB2;
    }
    text {
      font-size: 10px;
      font-family: sans-serif;
      /* in elk's coordinates "hanging" would be the correct value" */
      dominant-baseline: hanging;
      text-align: left;
    }
    g.port > text {
      font-size: 8px;
    }
    polyline {
      fill: none;
      stroke: black;
      stroke-width: 1;
    }
    path {
      fill: none;
      stroke: black;
      stroke-width: 1;
    }
  `;
  this._defs = new Xml(
    "marker",
    {
      "id": "arrow",
      "markerWidth": "10",
      "markerHeight": "8",
      "refX": "10",
      "refY": "4",
      "orient": "auto"
      },
      [
        new Xml("path", {"d": "M0,7 L10,4 L0,1 L0,7", "style": "fill: #000000;"})
      ]
  );

  this.reset();
}

Renderer.prototype = {
  constructor: Renderer,

  reset() {
    // internal housekeeping
    this._edgeRoutingStyle = {
      __global: "POLYLINE"
    };
    this._parentIds = {};
    this._edgeParents = {};
  },

  init(root) {
    // reset
    this.reset();
    this.registerParentIds(root);
    this.registerEdges(root);
  },

  /* Utility methods. */

  // edges can be specified anywhere, there coordinates however are relative
  //  a) to the source node's parent
  //  b) the source node, if the target is a descendant of the source node
  isDescendant(parent, node) {
    var current = node.id;
    while (this._parentIds[current]) {
      current = this._parentIds[current];
      if (current == parent.id) {
          return true;
      }
    }
    return false;
  },

  getOption(e, id) {
    if (!e) {
      return undefined;
    }
    if (e.id) {
      return e.id;
    }
    var suffix = id.substr(id.lastIndexOf('.') + 1);
    if (e[suffix]) {
      return e[suffix];
    }
    return undefined;
  },

  registerParentIds(p) {
    this._edgeParents[p.id] = [];
    if (p.properties) {
      var er = this.getOption(p.properties);
      if (er) {
        this._edgeRoutingStyle[p.id] = er;
      }
    }
    (p.children || []).forEach((c) => {
      this._parentIds[c.id] = p.id;
      this.registerParentIds(c);
    });
  },

  registerEdges(p) {
    (p.edges || []).forEach((e) => {
      e.sources.forEach(source_id => {
        e.targets.forEach(target_id => {
          if (source_id.includes(":")) {
            source_id = source_id.slice(0, source_id.indexOf(":"));
          }
          if (!this.isDescendant(source_id, target_id)) {
            source_id = this._parentIds[source_id];
          }
          this._edgeParents[source_id].push(e);
        });
      });
    });
    (p.children || []).forEach(c => this.registerEdges(c));
  },

  /*
   * Rendering methods.
   */

  renderRoot(root, styles="DEFAULT", defs="DEFAULT") {
    var children = [];

    var defsChildren = [];
    if (styles != null || defs != null) {
      if (styles != null) {
        defsChildren.push(this.svgCss(styles == "DEFAULT"? this._style: styles));
      }
      if (defs != null) {
        defsChildren.push(defs == "DEFAULT"? this._defs: new Text(defs));
      }
      children.push(new Xml("defs", {}, defsChildren));
    }

    children.push(this.renderGraph(root));

    return new Xml(
      "svg",
      {
          "version": "1.1",
          "xmlns": "http://www.w3.org/2000/svg",
          "width": root.width || 100,
          "height": root.height || 100
      },
      children
    );
  },

  renderGraph(graph) {
    var children = [];

    // paint edges first so that ports are drawn on top of them
    for (const edge of this._edgeParents[graph.id]) {
      children.push(this.renderEdge(edge, graph));
      if (edge.labels) {
        children.push(this.renderLabels(edge.labels));
      }
    }
    for (const child of graph.children) {
      children.push(this.renderRect(child));

      if (child.ports || child.labels) {
        children.push(this.renderPortsAndLabels(child))
      }

      if (child.children != null && child.children.length > 0) {
        children.push(this.renderGraph(child));
      }
    }
    var properties = {};
    if (graph.x || graph.y) {
      properties["transform"] = `translate(${graph.x || 0},${graph.y || 0})`;
    }
    return new Xml("g", properties, children);
  },

  renderPortsAndLabels(node) {
    var children = [];

    for (const p of node.ports) {
      children.push(this.renderRect(p));
      if (p.labels) {
        children.push(this.renderPort(p));
      }
    }
    if (node.labels) {
      for (const l of node.labels) {
        children.push(this.renderLabel(l));
      }
    }

    return new Xml("g", {"transform": `translate(${node.x || 0},${node.y || 0})`}, children);
  },

  renderRect(node) {
    return new Xml("rect", {
      ...this.idClass(node, "node"),
      ...this.posSize(node),
      ...this.style(node),
      ...this.attributes(node)
    })
  },

  renderPort(port) {
    return new Xml(
      "g",
      {
          ...this.idClass(port, "port"),
          "transform": `translate(${port.x || 0},${port.y || 0})`
      },
      this.renderLabels(port.labels)
    )
  },

  renderEdge(edge, node) {
    var bends = this.getBends(edge.sections);

    if (this._edgeRoutingStyle[node.id] == "SPLINES" || this._edgeRoutingStyle.__global == "SPLINES") {
      return this.renderPath(edge, bends);
    }
    return this.renderPolyline(edge, bends);
  },

  renderPath(edge, bends) {
    return new Xml("path", {
      "d": this.bendsToSpline(bends),
      ...this.idClass(edge, "edge"),
      ...this.style(edge),
      ...this.attributes(edge),
    });
  },

  renderPolyline(edge, bends) {
    return new Xml("polyline", {
      "points": this.bendsToPolyline(bends),
      ...this.idClass(edge, "edge"),
      ...this.style(edge),
      ...this.attributes(edge),
    });
  },

  getBends(sections) {
    var bends = [];
    if (sections && sections.length > 0) {
      sections.forEach(section => {
        if (section.startPoint) {
          bends.push(section.startPoint);
        }
        if (section.bendPoints) {
          bends = bends.concat(section.bendPoints);
        }
        if (section.endPoint) {
          bends.push(section.endPoint);
        }
      });
    }
    return bends;
  },

  renderLabels(labels) {
    return (labels || []).map(l => this.renderLabel(l))
  },

  renderLabel(label) {
    return new Xml("text", {
      ...this.idClass(label, "label"),
      ...this.posSize(label),
      ...this.style(label),
      ...this.attributes(label),
    }, [
      new Text(label.text)
    ]);
  },

  bendsToPolyline(bends) {
    return bends.map(bend => `${bend.x},${bend.y}`).join(" ")
  },

  bendsToSpline(bends) {
    if (!bends.length) {
      return ""
    }

    let {x, y} = bends[0];
    points = [`M${x} ${y}`]

    for (let i = 1; i < bends.length; i = i+3) {
      var left = bends.length - i;
      if (left == 1) {
        points.push(`L${bends[i].x + " " + bends[i].y}`);
      } else if (left == 2) {
        points.push(`Q${bends[i].x + " " + bends[i].y}`);
        points.push(bends[i+1].x + " " + bends[i+1].y);
      } else {
        points.push(`C${bends[i].x + " " + bends[i].y}`);
        points.push(bends[i+1].x + " " + bends[i+1].y);
        points.push(bends[i+2].x + " " + bends[i+2].y);
      }
    }
    return points.join(" ");
  },

  svgCss(css) {
    if (css == "") {
      return "";
    }
    return new Xml("style", {"type": "text/css"}, [
      new Cdata(
        new Text(css)
      )
    ]);
  },

  posSize(e) {
    return {
      "x": e.x || 0,
      "y": e.y || 0,
      "width": e.width || 0,
      "height": e.height || 0,
    };
  },

  idClass(e, className) {
    var elemClasses = Array.isArray(e.class)? e.class.join(" "): e.class;
    var classes = [elemClasses, className].filter(c => c).join(" ")

    var properties = {}
    if (e.id) {
      properties["id"] = e.id;
    }
    if (classes) {
      properties["class"] = classes;
    }
    return properties;
  },

  style(e) {
    if (!e.style) {
      return {}
    }
    return {
      "style": e.style
    }
  },

  attributes(e) {
    return e.attributes;
  },


  /*
   * Public API
   */

  toSvg(json, styles="DEFAULT", defs="DEFAULT") {
   this.init(json);
   var tree = this.renderRoot(json, styles, defs);
   return tree.render();
  }
};


window.elkSvg = {
  Renderer
};
