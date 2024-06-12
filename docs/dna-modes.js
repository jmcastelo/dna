importScripts("https://cdn.jsdelivr.net/pyodide/v0.25.0/full/pyodide.js");

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
  const env_spec = ['https://cdn.holoviz.org/panel/wheels/bokeh-3.4.1-py3-none-any.whl', 'https://cdn.holoviz.org/panel/1.4.4/dist/wheels/panel-1.4.4-py3-none-any.whl', 'pyodide-http==0.2.1', 'matplotlib', 'numpy', 'pandas']
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
  \nimport asyncio\n\nfrom panel.io.pyodide import init_doc, write_doc\n\ninit_doc()\n\nfrom numpy import zeros, array, linalg, pi, sqrt, abs, sin, cos, delete\nimport matplotlib.pyplot as plt\nfrom matplotlib.figure import Figure\nfrom matplotlib import colormaps\nimport pandas as pd\nimport panel as pn\nimport re\n\n\n\npn.extension(notifications = True)\nnotifications = pn.state.notifications\npn.extension('tabulator')\n\n\n\n# Data structures\n\nfreqs_xy, modes_xy = None, None\nfreqs_z, modes_z = None, None\n\n\n\n# Input and output widgets\n\ncard_margin = 20\n\ndna_sequence_input = pn.widgets.TextInput(name = 'DNA sequence', placeholder = 'DNA sequence...', sizing_mode = 'stretch_width')\ndna_sequence_length_indicator = pn.indicators.Number(name = 'Length', value = 0, sizing_mode = 'stretch_width')\na_length_indicator = pn.indicators.Number(name = '#A', value = 0, sizing_mode = 'stretch_width')\nt_length_indicator = pn.indicators.Number(name = '#T', value = 0, sizing_mode = 'stretch_width')\nc_length_indicator = pn.indicators.Number(name = '#C', value = 0, sizing_mode = 'stretch_width')\ng_length_indicator = pn.indicators.Number(name = '#G', value = 0, sizing_mode = 'stretch_width')\n\n\ndistance_input = pn.widgets.FloatInput(name = 'Distance (Angstrom)', value = 3.4, step = 0.1, start = 1.0e-6)\ntwist_angle_input = pn.widgets.FloatInput(name = 'Twist angle (rad)', value = pi / 5, step = 0.01, start = 0.0, end = 2 * pi)\n\nplot_width_input = pn.widgets.FloatInput(name = 'Plot width (inches)', value = 8, step = 0.1, start = 0.1)\nplot_height_input = pn.widgets.FloatInput(name = 'Plot height (inches)', value = 6, step = 0.1, start = 0.1)\n\nplot_button = pn.widgets.Button(name = 'Plot selected vibrational modes', button_type = 'primary', margin = card_margin)\n\n\n\n# Output tables\n\nfreqs_xy_table = pn.widgets.Tabulator(show_index = False, disabled = True, selectable = 'checkbox', pagination = 'local', page_size = 10, align = ('start'))\nfreqs_z_table = pn.widgets.Tabulator(show_index = False, disabled = True, selectable = 'checkbox', pagination = 'local', page_size = 10, align = ('start'))\n\n\n\n# Plot and HTML panes\n\npane_margin = 50\n\nfreqs_pane = pn.pane.Matplotlib(align = ('center', 'start'), sizing_mode = 'stretch_width', margin = pane_margin)\nmodes_pane = pn.pane.Matplotlib(align = ('center', 'start'), sizing_mode = 'stretch_width', margin = pane_margin)\n\nlinks_pane = pn.pane.HTML('''<a href="https://doi.org/10.1016/j.jtbi.2015.11.018">Reference</a>''')\n\n\n# Colormaps\n\ncmaps = dict()\ncmaps_list = list(colormaps)\nfor c in cmaps_list:\n    cmaps[c] = plt.get_cmap(c)\ncmaps_picker = pn.widgets.ColorMap(name = 'Colormap', options = cmaps, value_name = 'rainbow', swatch_width = 200)\n\n\n\n# Check letters\n\ndef check_dna_sequence(event):\n    # Check validity of letters\n\n    valid_letters = 'ATCGatcg'\n    input_dna_sequence = event.new\n\n    for i, letter in enumerate(input_dna_sequence):\n        if not letter in valid_letters:\n            notifications.error(f"Invalid base '{letter}' found at position {i + 1}", duration = 7000)\n\n\n    # Delete invalid letters\n    dna_sequence = re.sub('[^%s]' % valid_letters, '', input_dna_sequence)\n\n    # Count nucleotids\n\n    dna_sequence_lower = dna_sequence.lower();\n\n    num_a = dna_sequence_lower.count('a')\n    num_t = dna_sequence_lower.count('t')\n    num_c = dna_sequence_lower.count('c')\n    num_g = dna_sequence_lower.count('g')\n\n    # Update\n\n    event.obj.value_input = dna_sequence_lower\n\n    dna_sequence_length_indicator.value = len(dna_sequence_lower)\n\n    a_length_indicator.value = num_a\n    t_length_indicator.value = num_t\n    c_length_indicator.value = num_c\n    g_length_indicator.value = num_g\n\n\n\ndna_sequence_input.param.watch(check_dna_sequence, 'value_input')\n\n\n\n# Construct dynamic matrices\n\ndef set_dynamic_matrices():\n    # Main variables\n    dna_sequence = dna_sequence_input.value.upper()\n    distance = distance_input.value\n    twist_angle = twist_angle_input.value\n\n    # Conversion factors\n    ang_to_meter = 1.0e-10\n\n    # Fundamental constants\n    qe = 1.602176634e-19\n    me = 9.1093837015e-31\n    e0 = 8.8541878128e-12\n\n    # DNA base pair electronic frequencies\n    omega = dict()\n\n    omega['A', 'xx'] = 3.062e15\n    omega['A', 'yy'] = 2.822e15\n    omega['A', 'zz'] = 4.242e15\n\n    omega['T', 'xx'], omega['T', 'yy'], omega['T', 'zz'] = omega['A', 'xx'], omega['A', 'yy'], omega['A', 'zz']\n\n    omega['C', 'xx'] = 3.027e15\n    omega['C', 'yy'] = 2.722e15\n    omega['C', 'zz'] = 4.244e15\n\n    omega['G', 'xx'], omega['G', 'yy'], omega['G', 'zz'] = omega['C', 'xx'], omega['C', 'yy'], omega['C', 'zz']\n\n\n    # Potential matrix elements\n\n    gamma = dict()\n\n    gamma['xx'] = (cos(twist_angle) * qe ** 2) / (4 * pi * e0 * (distance * ang_to_meter) ** 3)\n    gamma['yy'] = gamma['xx']\n    gamma['zz'] = -(qe ** 2) / (2 * pi * e0 * (distance * ang_to_meter) ** 3)\n    gamma['xy'] = -(sin(twist_angle) * qe ** 2) / (4 * pi * e0 * (distance * ang_to_meter) ** 3)\n\n    # DNA Sequence length\n    sequence_length = len(dna_sequence)\n\n    # Longitudinal potential matrix\n\n    Vz = zeros([sequence_length, sequence_length])\n\n    for i in range(sequence_length):\n        Vz[i, i] = me * omega[dna_sequence[i], 'zz'] ** 2\n\n    for i in range(sequence_length - 1):\n        Vz[i, i + 1] = gamma['zz']\n        Vz[i + 1, i] = gamma['zz']\n\n    # Transverse potential matrix\n\n    Vxy = zeros([2 * sequence_length, 2 * sequence_length])\n\n    for i in range(sequence_length):\n        Vxy[i, i] = me * omega[dna_sequence[i], 'xx'] ** 2\n        Vxy[sequence_length + i, sequence_length + i] = me * omega[dna_sequence[i], 'yy'] ** 2\n\n    for i in range(sequence_length - 1):\n        Vxy[i, i + 1] = gamma['xx']\n        Vxy[i + 1, i] = gamma['xx']\n        Vxy[sequence_length + i, sequence_length + i + 1] = gamma['yy']\n        Vxy[sequence_length + i + 1, sequence_length + i] = gamma['yy']\n\n    for i in range(sequence_length - 1):\n        Vxy[i, sequence_length + i + 1] = gamma['xy']\n        Vxy[sequence_length + i + 1, i] = gamma['xy']\n        Vxy[i + 1, sequence_length + i] = -gamma['xy']\n        Vxy[sequence_length + i, i + 1] = -gamma['xy']\n\n    return Vxy, Vz\n\n\n\n# Process egivenvalues and eigenvectors and get relevant quantities\n\ndef get_frequencies_and_modes(eigenvalues, eigenvectors, transverse):\n    # Fundamental constants\n\n    ev = 1.602176634e-19 # Electronvolt\n    hbar = 6.62607015e-34 / (2 * pi) # Planck's constant\n    e_PO = 0.23 # Reference Energy\n    me = 9.1093837015e-31 # Electron mass\n\n    # Keep positive eigenvalues only\n\n    negative = []\n    for i in range(len(eigenvalues)):\n        if eigenvalues[i] < 0:\n            negative.append(i)\n\n    freqs = delete(eigenvalues, negative)\n    modes = delete(eigenvectors, negative, 1)\n\n    # Compute quantities of interest\n\n    frequencies = sqrt(freqs / me) * hbar / (2.0 * e_PO * ev)\n\n    if transverse:\n        half_comp = int(modes.shape[0] / 2)\n        mode_amplitudes = array([[abs(modes[comp, vec]) ** 2 + abs(modes[half_comp + comp, vec]) ** 2 for comp in range(half_comp)] for vec in range(modes.shape[1])])\n    else:\n        mode_amplitudes = array([[abs(modes[comp, vec]) ** 2 for comp in range(modes.shape[0])] for vec in range(modes.shape[1])])\n\n    return frequencies, mode_amplitudes\n\n\n\n# Plot frequencies/energies\n\ndef plot_frequencies():\n    # Construct tags and indexes for selected transverse and longitudinal modes\n\n    sel_idx_xy = freqs_xy_table.selection\n    sel_idx_z = freqs_z_table.selection\n\n    if len(sel_idx_xy) + len(sel_idx_z) <= 0:\n        return\n\n    sel_idx_xy.sort()\n    sel_idx_z.sort()\n\n    tag_xy = [f"T{i}" for i in sel_idx_xy]\n    tag_z = [f"L{i}" for i in sel_idx_z]\n\n    tag = tag_xy + tag_z\n\n    # Select frequencies\n\n    sel_freqs_xy = freqs_xy[sel_idx_xy]\n    sel_freqs_z = freqs_z[sel_idx_z]\n\n    # Construct indexes for colormap\n\n    max_len = max(len(freqs_xy) - 1, len(freqs_z) - 1)\n    norm_idx_xy = [s / max_len for s in sel_idx_xy]\n    norm_idx_z = [s / max_len for s in sel_idx_z]\n    norm_idx = norm_idx_xy + norm_idx_z\n\n    # Join lists\n    freqs = sel_freqs_xy.tolist() + sel_freqs_z.tolist()\n\n    # Plot options\n\n    fig = Figure(figsize = (plot_width_input.value, plot_height_input.value), layout = 'constrained')\n    ax = fig.subplots()\n\n    ax.set_xlabel('Modes')\n    ax.set_ylabel('Energy (Epo)')\n\n    color_map = plt.get_cmap(cmaps_picker.value_name)\n\n    ax.bar(tag, freqs, color = color_map(norm_idx))\n\n    # Remove tick labels to enhance visibility\n\n    tick_labels = ax.xaxis.get_majorticklabels()\n    if len(tick_labels) >= 20:\n        for i, label in enumerate(tick_labels):\n            if i % int(len(tick_labels) / 20) > 0:\n                label.set_visible(False)\n\n    return fig\n\n\n\n# Plot modes\n\ndef plot_modes():\n    # Construct tags and markers for selected transverse and longitudinal modes\n\n    sel_idx_xy = freqs_xy_table.selection\n    sel_idx_z = freqs_z_table.selection\n\n    if len(sel_idx_xy) + len(sel_idx_z) <= 0:\n        return\n\n    sel_idx_xy.sort()\n    sel_idx_z.sort()\n\n    tag_xy = [f"T{i}" for i in sel_idx_xy]\n    tag_z = [f"L{i}" for i in sel_idx_z]\n\n    tag = tag_xy + tag_z\n\n    marker_xy = ['.-' for i in sel_idx_xy]\n    marker_z = ['.:' for i in sel_idx_z]\n    marker = marker_xy + marker_z\n\n    # Select modes\n\n    sel_modes_xy = modes_xy[sel_idx_xy]\n    sel_modes_z = modes_z[sel_idx_z]\n\n    # Construct indexes for colormap\n\n    max_len = max(len(modes_xy) - 1, len(modes_z) - 1)\n    norm_idx_xy = [s / max_len for s in sel_idx_xy]\n    norm_idx_z = [s / max_len for s in sel_idx_z]\n    norm_seq = norm_idx_xy + norm_idx_z\n\n    # Join lists\n    modes = sel_modes_xy.tolist() + sel_modes_z.tolist()\n\n    # Plot options\n\n    fig = Figure(figsize = (plot_width_input.value, plot_height_input.value), layout = 'constrained')\n    ax = fig.subplots()\n\n    ax.set_xlabel('DNA Sequence')\n    ax.set_ylabel('Amplitude')\n\n    color_map = plt.get_cmap(cmaps_picker.value_name)\n\n    for i, mode in enumerate(modes):\n        ax.plot(mode, marker[i], color = color_map(norm_seq[i]), label = tag[i])\n\n    dna_sequence = dna_sequence_input.value.upper()\n    nucleotids = [letter for letter in dna_sequence]\n    ax.set_xticks(range(len(nucleotids)), nucleotids)\n\n    ax.legend(loc = 'upper left', bbox_to_anchor = (0.05, -0.2), ncol = 10, fontsize = 'x-small')\n\n    return fig\n\n\n\n# Frequencies/energies tables\n\ndef set_freqs_xy_table_data():\n    # Tags\n\n    idx = list(range(len(freqs_xy)))\n    tag = [f"T{i:03d}" for i in idx]\n\n    # Data\n\n    df = pd.DataFrame({\n        'Modes': tag,\n        'Energies (Epo)': freqs_xy.tolist()\n    }, index = idx)\n\n    freqs_xy_table.value = df\n\n\n\ndef set_freqs_z_table_data():\n    # Tags\n\n    idx = list(range(len(freqs_z)))\n    tag = [f"L{i:03d}" for i in idx]\n\n    # Data\n\n    df = pd.DataFrame({\n        'Modes': tag,\n        'Energies (Epo)': freqs_z.tolist()\n    }, index = idx)\n\n    freqs_z_table.value = df\n\n\n\n# Set plots\n\ndef set_outputs(event):\n    fig_freqs = plot_frequencies()\n    freqs_pane.object = fig_freqs\n\n    fig_modes = plot_modes()\n    modes_pane.object = fig_modes\n\n\n\npn.bind(set_outputs, plot_button, watch = True)\n\n\n\n# Remove plots\n\ndef clear_outputs():\n    freqs_pane.object = None\n    modes_pane.object = None\n\n\n\n# Main computing function\n\ndef compute_vibrational_modes(event):\n    dna_sequence = dna_sequence_input.value.upper()\n\n    sequence_length = len(dna_sequence)\n\n    # Checks\n\n    if sequence_length < 2:\n        notifications.error('Length of DNA chain smaller than 2!')\n        return\n\n    # Build dynamic matrices\n\n    Vxy, Vz = set_dynamic_matrices()\n\n    # Diagonalize dynamic matrices\n\n    eigenvalues_xy, eigenvectors_xy = linalg.eigh(Vxy)\n    eigenvalues_z, eigenvectors_z = linalg.eigh(Vz)\n\n    # Get frequencies and modes\n\n    global freqs_xy, modes_xy\n    global freqs_z, modes_z\n\n    freqs_xy, modes_xy = get_frequencies_and_modes(eigenvalues_xy, eigenvectors_xy, True)\n    freqs_z, modes_z = get_frequencies_and_modes(eigenvalues_z, eigenvectors_z, False)\n\n    set_freqs_xy_table_data()\n    set_freqs_z_table_data()\n\n    # Clear plots\n    clear_outputs()\n\n\n\ndna_sequence_input.param.watch(compute_vibrational_modes, 'value')\ndistance_input.param.watch(compute_vibrational_modes, 'value')\ntwist_angle_input.param.watch(compute_vibrational_modes, 'value')\n\n\n\n# Output layouts\n\ncounts_layout = pn.Row(dna_sequence_length_indicator, a_length_indicator, t_length_indicator, c_length_indicator, g_length_indicator, sizing_mode = 'stretch_width')\npanes_layout = pn.Column(freqs_pane, modes_pane, sizing_mode = 'stretch_both')\n\n# Controls layout\n\ncomputation_parameters_card = pn.Card(distance_input, twist_angle_input, title = 'Parameters', margin = card_margin)\nxy_modes_card = pn.Card(freqs_xy_table, title = 'Transverse modes', margin = card_margin)\nt_modes_card = pn.Card(freqs_z_table, title = 'Longitudinal modes', margin = card_margin)\nplot_options_card = pn.Card(plot_width_input, plot_height_input, cmaps_picker, title = 'Plot options', margin = card_margin)\n\ncontrols_layout = pn.FlexBox(xy_modes_card, t_modes_card, computation_parameters_card, plot_options_card, plot_button, flex_direction = 'row', sizing_mode = 'stretch_width')\n\n# Main layout\n\nmain_layout = pn.Column(dna_sequence_input, counts_layout,  pn.layout.Divider(), controls_layout, pn.layout.Divider(), panes_layout, sizing_mode = 'stretch_both', scroll = True)\n\n# Template\n\ntemplate = pn.template.BootstrapTemplate(title = 'DNA Electronic Quantum Vibrational Modes', header = links_pane)\ntemplate.main.append(main_layout)\ntemplate.servable()\n\n\nawait write_doc()
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
    from panel.io.pyodide import _convert_json_patch
    state.curdoc.apply_json_patch(_convert_json_patch(patch), setter='js')
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