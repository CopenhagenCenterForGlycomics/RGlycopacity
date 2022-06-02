
if (!Array.prototype.flat) {
  Object.defineProperty(Array.prototype, 'flat', {
    configurable: true,
    value: function flat () {
      var depth = isNaN(arguments[0]) ? 1 : Number(arguments[0]);

      return depth ? Array.prototype.reduce.call(this, function (acc, cur) {
        if (Array.isArray(cur)) {
          acc.push.apply(acc, flat.call(cur, depth - 1));
        } else {
          acc.push(cur);
        }

        return acc;
      }, []) : Array.prototype.slice.call(this);
    },
    writable: true
  });
}

const set_colours = async function(topiary,colours) {
  for (let exp_targ of topiary.querySelectorAll('g[id*="HGNC"]')) {
    let exp_id = exp_targ.id;
    let [_,hgnc,entrez] = exp_id.match(/HGNC:([^_]+)_entrez:(\d+)/);
    entrez = parseInt(entrez);
    let expression_row;
    if (colours[hgnc]) {
      expression_row = { colour: colours[hgnc]};
    }
    if ( ! expression_row ) {
      continue;
    }
    let target_colour = expression_row.colour || 'inherit';
    if (exp_targ.querySelector('path')) {
      if (target_colour === 'inherit') {
        delete exp_targ.querySelector('path').style.fill;
        delete exp_targ.querySelector('path').style.opacity;
        if (exp_targ.querySelector('path').style.removeProperty) {
          exp_targ.querySelector('path').style.removeProperty('fill');
          exp_targ.querySelector('path').style.removeProperty('opacity');
        }
      } else {
        exp_targ.querySelector('path').style.fill = `#${target_colour.replace('#','')}`;
        exp_targ.querySelector('path').style.opacity = 1;
      }
    }
    if (exp_targ.querySelector('circle')) {
      if (target_colour === 'inherit') {
        delete exp_targ.querySelector('circle').style.fill;
        delete exp_targ.querySelector('circle').style.opacity;
        if (exp_targ.querySelector('circle').style.removeProperty) {
          exp_targ.querySelector('circle').style.removeProperty('fill');
          exp_targ.querySelector('circle').style.removeProperty('opacity');
        }
      } else {
        exp_targ.querySelector('circle').style.fill = `#${target_colour.replace('#','')}`;
        exp_targ.querySelector('circle').style.opacity = 1;
      }
    }
  }
}

const set_size = async function(svg_doc,portrait) {
  let renderer = load_sugar_svg(svg_doc);
  renderer.element.canvas.firstElementChild.setAttribute('width',portrait ? '594mm' : '840mm');
  renderer.element.canvas.firstElementChild.setAttribute('height',portrait ? '840mm' : '594mm');
  console.r.assign('svg_temp',renderer.element.canvas.innerHTML);
}

const recolour_labels = async function(svg_doc,colours) {
  let renderer = load_sugar_svg(svg_doc);
  await set_colours(renderer.element.canvas,colours);
  renderer.element.canvas.firstElementChild.setAttribute('width','594mm');
  renderer.element.canvas.firstElementChild.setAttribute('height','840mm');
  console.r.assign('svg_temp',renderer.element.canvas.innerHTML)
}

const filter_entries = function(object,keys) {
  let wanted_keys = Object.keys(object).filter( key => keys.indexOf(key) >= 0);
  let result = {};
  for (let key of wanted_keys) {
    result[key] = object[key];
  }
  return result;
}

const perform_glycotopiary = async function(svg_doc,reactions_json,wanted_genes) {
  try {
    let reactions = JSON.parse(reactions_json);
    if (wanted_genes === null || !Array.isArray(wanted_genes) )  {
      wanted_genes = Object.keys(reactions); 
    }
    reactions = filter_entries(reactions,wanted_genes);
    let renderer = load_sugar_svg(svg_doc) 
    tag_supported(renderer,reactions);
    let wanted_symbol = renderer.groupTag;
    for (let sugar of renderer.sugars) {
      let unsupported = [...sugar.composition()].filter( res => !res.getTag(wanted_symbol) );
      for (let res of unsupported) {
        if ( ! renderer.rendered.get(res)) {
          continue;
        }
        let { residue, linkage } = renderer.rendered.get(res);
        residue.element.style.opacity = 0.1;
        if (linkage) {
          linkage.style.opacity = 0.1;
        }
      }
    }
    console.r.assign('svg_temp',renderer.element.canvas.innerHTML)
  } catch (err) {
    console.log(err.stack);
  }
}

const lookup_opacity = function(root, id) {
  id = id.replace('url(','').replace(')','');
  let mask = root.querySelector(id);
  return mask.querySelector('rect').style.fillOpacity;
}

const update_opacity = function(el, opacity) {
  for (let child of el.querySelectorAll('path')) {
    child.style.opacity = opacity;
  }
}

const fix_opacities = async function(svg_doc) {
  let renderer = load_sugar_svg(svg_doc);
  let root = renderer.element.canvas;
  for (let mask_el of root.querySelectorAll('use[mask]')) {
    let opacity = lookup_opacity(root,mask_el.getAttribute('mask'));
    update_opacity(root.querySelector(mask_el.getAttribute('xlink:href')),opacity);
  }
  console.r.assign('svg_temp',renderer.element.canvas.innerHTML);
}