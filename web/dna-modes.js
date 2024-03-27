importScripts("https://cdn.jsdelivr.net/pyodide/v0.24.1/full/pyodide.js");

function sendPatch(patch, buffers, msg_id) {
  self.postMessage({
    type: 'patch',
    patch: patch,
    buffers: buffers
  })
}

async function startApplication() {
  console.log("Loading pyodide!");
  self.postMessage({type: 'status', msg: 'Loading pyodide'})
  self.pyodide = await loadPyodide();
  self.pyodide.globals.set("sendPatch", sendPatch);
  console.log("Loaded!");
  await self.pyodide.loadPackage("micropip");
  const env_spec = ['https://cdn.holoviz.org/panel/wheels/bokeh-3.3.4-py3-none-any.whl', 'https://cdn.holoviz.org/panel/1.3.8/dist/wheels/panel-1.3.8-py3-none-any.whl', 'pyodide-http==0.2.1', 'matplotlib', 'numpy', 'pandas']
  for (const pkg of env_spec) {
    let pkg_name;
    if (pkg.endsWith('.whl')) {
      pkg_name = pkg.split('/').slice(-1)[0].split('-')[0]
    } else {
      pkg_name = pkg
    }
    self.postMessage({type: 'status', msg: `Installing ${pkg_name}`})
    try {
      await self.pyodide.runPythonAsync(`
        import micropip
        await micropip.install('${pkg}');
      `);
    } catch(e) {
      console.log(e)
      self.postMessage({
	type: 'status',
	msg: `Error while installing ${pkg_name}`
      });
    }
  }
  console.log("Packages loaded!");
  self.postMessage({type: 'status', msg: 'Executing code'})
  const code = `
  
import asyncio

from panel.io.pyodide import init_doc, write_doc

init_doc()

from numpy import array, zeros, linalg, arange, pi, arange, loadtxt, sqrt, dot, abs, sin, cos, delete, concatenate
import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.cm import get_cmap
from math import isnan
import pandas as pd
import panel as pn


pn.extension(notifications = True)
pn.extension('tabulator')



def set_dynamic_matrices(dna_chain, distance = 3.4, twist_angle = pi / 5):
    # Conversion factors
    rad_to_Hz = 1.0 #/ (2 * pi)
    ang_to_meter = 1.0e-10

    # DNA base pair electronic frequencies
    omega = dict()

    omega['A', 'xx'] = 3.062e15 * rad_to_Hz;
    omega['A', 'yy'] = 2.822e15 * rad_to_Hz;
    omega['A', 'zz'] = 4.242e15 * rad_to_Hz;

    omega['T', 'xx'], omega['T', 'yy'], omega['T', 'zz'] = omega['A', 'xx'], omega['A', 'yy'], omega['A', 'zz']

    omega['C', 'xx'] = 3.027e15 * rad_to_Hz;
    omega['C', 'yy'] = 2.722e15 * rad_to_Hz;
    omega['C', 'zz'] = 4.244e15 * rad_to_Hz;

    omega['G', 'xx'], omega['G', 'yy'], omega['G', 'zz'] = omega['C', 'xx'], omega['C', 'yy'], omega['C', 'zz']

    # Fundamental constants
    qe = 1.602176634e-19
    me = 9.1093837015e-31
    e0 = 8.8541878128e-12

    # Potential matrix elements
    gamma = dict()

    gamma['xx'] = (cos(twist_angle) * qe ** 2) / (4 * pi * e0 * me * (distance * ang_to_meter) ** 3)
    gamma['yy'] = gamma['xx']
    gamma['zz'] = -(qe ** 2) / (2 * pi * e0 * me * (distance * ang_to_meter) ** 3)
    gamma['xy'] = -(sin(twist_angle) * qe ** 2) / (4 * pi * e0 * me * (distance * ang_to_meter) ** 3)
    gamma['yx'] = -gamma['xy']

    # DNA Chain length
    n_chain = len(dna_chain)

    # Longitudinal potential matrix
    Vz = zeros([n_chain, n_chain])

    for i in range(n_chain):
        Vz[i, i] = omega[dna_chain[i], 'zz'] ** 2

    for i in range(n_chain - 1):
        Vz[i, i + 1] = gamma['zz']
        Vz[i + 1, i] = gamma['zz']

    # Transverse potential matrix
    Vxy = zeros([2 * n_chain, 2 * n_chain])

    for i in range(n_chain):
        Vxy[i, i] = omega[dna_chain[i], 'xx'] ** 2
        Vxy[n_chain + i, n_chain + i] = omega[dna_chain[i], 'yy'] ** 2

    for i in range(n_chain - 1):
        Vxy[i, i + 1] = gamma['xx']
        Vxy[i + 1, i] = gamma['xx']
        Vxy[n_chain + i, n_chain + i + 1] = gamma['yy']
        Vxy[n_chain + i + 1, n_chain + i] = gamma['yy']

    for i in range(n_chain - 1):
        Vxy[i, n_chain + i + 1] = gamma['xy']
        Vxy[n_chain + i + 1, i] = gamma['xy']
        Vxy[i + 1, n_chain + i] = -gamma['xy']
        Vxy[n_chain + i, i + 1] = -gamma['xy']

    return Vxy, Vz



def get_frequencies_and_modes(eigenvalues, eigenvectors):
    ev = 1.602176634e-19
    hbar = 6.62607015e-34 / (2 * pi)
    e_PO = 0.23 # Reference Energy

    negative = []
    for i in range(len(eigenvalues)):
        if eigenvalues[i] < 0:
            negative.append(i)

    freqs = delete(eigenvalues, negative)
    modes = delete(eigenvectors, negative, 0)

    frequencies = sqrt(freqs) * hbar / (2.0 * e_PO * ev)

    return frequencies, modes



# Colormaps
cmaps = dict()
cmaps_list = list(colormaps)
for c in cmaps_list:
    cmaps[c] = get_cmap(c)
cmaps_picker = pn.widgets.ColorMap(options = cmaps, value_name = 'rainbow')



def plot_frequencies(freqs_xy, freqs_z, dna_chain):
    seq_xy = list(range(len(freqs_xy)))
    seq_z = list(range(len(freqs_z)))

    tag_xy = [f"T{i}" for i in seq_xy]
    tag_z = [f"L{i}" for i in seq_z]

    tag = tag_xy + tag_z

    max_len = max(len(freqs_xy) - 1, len(freqs_z) - 1)
    norm_seq_xy = [s / max_len for s in seq_xy]
    norm_seq_z = [s / max_len for s in seq_z]
    norm_seq = norm_seq_xy + norm_seq_z

    freqs = freqs_xy.tolist() + freqs_z.tolist()

    # Plot options

    fig, ax = plt.subplots(figsize = (4, 4))

    ax.set_xlabel(dna_chain)

    color_map = plt.get_cmap(cmaps_picker.value_name)

    b = ax.bar(tag, freqs, color = color_map(norm_seq))
    ax.bar_label(b, label_type = 'edge', fontsize = 6, fmt = '{:,.2f}')

    return fig



freqs_table = pn.widgets.Tabulator(show_index = False, disabled = True, align = ('center'))

def set_freqs_table_data(freqs_xy, freqs_z):
    idx_xy = list(range(len(freqs_xy)))
    idx_z = list(range(len(freqs_z)))

    tag_xy = [f"T{i}" for i in idx_xy]
    tag_z = [f"L{i}" for i in idx_z]

    idx = idx_xy + idx_z
    tag = tag_xy + tag_z

    freqs = freqs_xy.tolist() + freqs_z.tolist()

    df = pd.DataFrame({
        'Modes': tag,
        'Energies (Epo)': freqs
    }, index = idx)

    freqs_table.value = df



# Input widgets
chain_input = pn.widgets.TextInput(name = 'DNA sequence:', placeholder = 'DNA sequence...')
distance_input = pn.widgets.FloatInput(name = 'Distance (Angstrom):', value = 3.4, step = 0.1, start = 1.0e-6)
twist_angle_input = pn.widgets.FloatInput(name = 'Twist angle (rad):', value = pi / 5, step = 0.01, start = 0.0, end = 2 * pi)
max_xy_mode_input = pn.widgets.IntInput(name = 'Max. transverse mode:', value = 2, start = 1)
max_z_mode_input = pn.widgets.IntInput(name = 'Max. longitudinal mode:', value = 2, start = 1)

freqs_pane = pn.pane.Matplotlib(align = ('center'))



def compute_vibrational_modes(event):
    dna_chain = chain_input.value.upper()

    n_chain = len(dna_chain)

    if n_chain < 2:
        pn.state.notifications.error('Length of DNA chain smaller than 2!')
        return

    chain_elements = 'ATCG'
    if not all(c in chain_elements for c in dna_chain):
        pn.state.notifications.error('Invalid chain: must be combination of A, T, C, G')
        return

    distance = distance_input.value
    twist_angle = twist_angle_input.value

    Vxy, Vz = set_dynamic_matrices(dna_chain, distance, twist_angle)

    eigenvalues_xy, eigenvectors_xy = linalg.eigh(Vxy)
    eigenvalues_z, eigenvectors_z = linalg.eigh(Vz)

    freqs_xy, modes_xy = get_frequencies_and_modes(eigenvalues_xy, eigenvectors_xy)
    freqs_z, modes_z = get_frequencies_and_modes(eigenvalues_z, eigenvalues_z)

    max_xy_mode = min(max_xy_mode_input.value, len(modes_xy))
    max_z_mode = min(max_z_mode_input.value, len(modes_z))

    max_xy_mode_input.value = max_xy_mode
    max_z_mode_input.value = max_z_mode

    freqs, modes = [], []
    for i in range(max_xy_mode):
        freqs.append(freqs_xy[i])
        modes.append(modes_xy[i])
    for i in range(max_z_mode):
        freqs.append(freqs_z[i])
        modes.append(modes_z[i])

    set_freqs_table_data(freqs_xy[0 : max_xy_mode], freqs_z[0 : max_z_mode])

    fig_freqs = plot_frequencies(freqs_xy[0 : max_xy_mode], freqs_z[0 : max_z_mode], dna_chain)
    freqs_pane.object = fig_freqs



compute_button = pn.widgets.Button(name = 'Compute vibrational modes', button_type = 'primary')
pn.bind(compute_vibrational_modes, compute_button, watch = True)

output_box = pn.FlexBox(freqs_table, freqs_pane, align_content = 'center', align_items = 'center', justify_content = 'space-evenly')

links_pane = pn.pane.HTML('''<a href="https://doi.org/10.1016/j.jtbi.2015.11.018">Reference</a>''')

cmaps_text = pn.widgets.StaticText(name = '', value = 'Colormap:')
cmaps_column = pn.Column(cmaps_text, cmaps_picker)

template = pn.template.BootstrapTemplate(title = 'DNA Electronic Quantum Vibrational Modes', sidebar = [chain_input, distance_input, twist_angle_input, max_xy_mode_input, max_z_mode_input, compute_button, cmaps_column, links_pane])
template.main.append(output_box)
template.servable()


await write_doc()
  `

  try {
    const [docs_json, render_items, root_ids] = await self.pyodide.runPythonAsync(code)
    self.postMessage({
      type: 'render',
      docs_json: docs_json,
      render_items: render_items,
      root_ids: root_ids
    })
  } catch(e) {
    const traceback = `${e}`
    const tblines = traceback.split('\n')
    self.postMessage({
      type: 'status',
      msg: tblines[tblines.length-2]
    });
    throw e
  }
}

self.onmessage = async (event) => {
  const msg = event.data
  if (msg.type === 'rendered') {
    self.pyodide.runPythonAsync(`
    from panel.io.state import state
    from panel.io.pyodide import _link_docs_worker

    _link_docs_worker(state.curdoc, sendPatch, setter='js')
    `)
  } else if (msg.type === 'patch') {
    self.pyodide.globals.set('patch', msg.patch)
    self.pyodide.runPythonAsync(`
    state.curdoc.apply_json_patch(patch.to_py(), setter='js')
    `)
    self.postMessage({type: 'idle'})
  } else if (msg.type === 'location') {
    self.pyodide.globals.set('location', msg.location)
    self.pyodide.runPythonAsync(`
    import json
    from panel.io.state import state
    from panel.util import edit_readonly
    if state.location:
        loc_data = json.loads(location)
        with edit_readonly(state.location):
            state.location.param.update({
                k: v for k, v in loc_data.items() if k in state.location.param
            })
    `)
  }
}

startApplication()