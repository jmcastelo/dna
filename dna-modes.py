from numpy import zeros, array, linalg, pi, sqrt, abs, sin, cos, delete
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib import colormaps
import pandas as pd
import panel as pn
import re
from io import BytesIO
import zipfile
from datetime import datetime



pn.extension(notifications = True)
notifications = pn.state.notifications
pn.extension('tabulator')



# Global arrays

Vxy, Vz = None, None

eigenvalues_xy, eigenvectors_xy = None, None
eigenvalues_z, eigenvectors_z = None, None

energies_xy, modes_xy = None, None
energies_z, modes_z = None, None



# Input and output widgets

dna_sequence_input = pn.widgets.TextInput(name = 'DNA sequence', placeholder = 'DNA sequence...', sizing_mode = 'stretch_width')
dna_sequence_length_indicator = pn.indicators.Number(name = 'Length', value = 0, sizing_mode = 'stretch_width')
a_length_indicator = pn.indicators.Number(name = '#A', value = 0, sizing_mode = 'stretch_width')
t_length_indicator = pn.indicators.Number(name = '#T', value = 0, sizing_mode = 'stretch_width')
c_length_indicator = pn.indicators.Number(name = '#C', value = 0, sizing_mode = 'stretch_width')
g_length_indicator = pn.indicators.Number(name = '#G', value = 0, sizing_mode = 'stretch_width')


distance_input = pn.widgets.FloatInput(name = 'Distance (Angstrom)', value = 3.4, step = 0.1, start = 1.0e-6)
twist_angle_input = pn.widgets.FloatInput(name = 'Twist angle (rad)', value = pi / 5, step = 0.01, start = 0.0, end = 2 * pi)

plot_width_input = pn.widgets.FloatInput(name = 'Plot width (inches)', value = 8, step = 0.1, start = 0.1)
plot_height_input = pn.widgets.FloatInput(name = 'Plot height (inches)', value = 6, step = 0.1, start = 0.1)

plot_button = pn.widgets.Button(name = 'Plot selected vibrational modes', disabled = True,  button_type = 'primary')
data_download_button = pn.widgets.FileDownload(label = 'Download input and output data', disabled = True, button_type = 'primary')



# Output tables

energies_xy_table = pn.widgets.Tabulator(show_index = False, disabled = True, selectable = 'checkbox', pagination = 'local', page_size = 10, align = ('start'))
energies_z_table = pn.widgets.Tabulator(show_index = False, disabled = True, selectable = 'checkbox', pagination = 'local', page_size = 10, align = ('start'))



# Plot and HTML panes

pane_margin = 50

energies_pane = pn.pane.Matplotlib(align = ('center', 'start'), sizing_mode = 'stretch_width', margin = pane_margin)
modes_pane = pn.pane.Matplotlib(align = ('center', 'start'), sizing_mode = 'stretch_width', margin = pane_margin)

links_pane = pn.pane.HTML('''<a href="https://doi.org/10.1016/j.jtbi.2015.11.018">Reference</a>''')


# Colormaps

cmaps = dict()
cmaps_list = list(colormaps)
for c in cmaps_list:
    cmaps[c] = plt.get_cmap(c)
cmaps_picker = pn.widgets.ColorMap(name = 'Colormap', options = cmaps, value_name = 'rainbow', swatch_width = 200)



# Check letters

def check_dna_sequence(event):
    # Check validity of letters

    valid_letters = 'ATCGatcg'
    input_dna_sequence = event.new

    for i, letter in enumerate(input_dna_sequence):
        if not letter in valid_letters:
            notifications.error(f"Invalid base '{letter}' found at position {i + 1}", duration = 7000)


    # Delete invalid letters
    dna_sequence = re.sub('[^%s]' % valid_letters, '', input_dna_sequence)

    # Count nucleotids

    dna_sequence_lower = dna_sequence.lower();

    num_a = dna_sequence_lower.count('a')
    num_t = dna_sequence_lower.count('t')
    num_c = dna_sequence_lower.count('c')
    num_g = dna_sequence_lower.count('g')

    # Update

    event.obj.value_input = dna_sequence_lower

    dna_sequence_length_indicator.value = len(dna_sequence_lower)

    a_length_indicator.value = num_a
    t_length_indicator.value = num_t
    c_length_indicator.value = num_c
    g_length_indicator.value = num_g



dna_sequence_input.param.watch(check_dna_sequence, 'value_input')



# Construct dynamic matrices

def set_dynamic_matrices():
    # Main variables
    dna_sequence = dna_sequence_input.value.upper()
    distance = distance_input.value
    twist_angle = twist_angle_input.value

    # Conversion factors
    ang_to_meter = 1.0e-10

    # Fundamental constants
    qe = 1.602176634e-19
    me = 9.1093837015e-31
    e0 = 8.8541878128e-12

    # DNA base pair electronic frequencies
    omega = dict()

    omega['A', 'xx'] = 3.062e15
    omega['A', 'yy'] = 2.822e15
    omega['A', 'zz'] = 4.242e15

    omega['T', 'xx'], omega['T', 'yy'], omega['T', 'zz'] = omega['A', 'xx'], omega['A', 'yy'], omega['A', 'zz']

    omega['C', 'xx'] = 3.027e15
    omega['C', 'yy'] = 2.722e15
    omega['C', 'zz'] = 4.244e15

    omega['G', 'xx'], omega['G', 'yy'], omega['G', 'zz'] = omega['C', 'xx'], omega['C', 'yy'], omega['C', 'zz']


    # Potential matrix elements

    gamma = dict()

    gamma['xx'] = (cos(twist_angle) * qe ** 2) / (4 * pi * e0 * (distance * ang_to_meter) ** 3)
    gamma['yy'] = gamma['xx']
    gamma['zz'] = -(qe ** 2) / (2 * pi * e0 * (distance * ang_to_meter) ** 3)
    gamma['xy'] = -(sin(twist_angle) * qe ** 2) / (4 * pi * e0 * (distance * ang_to_meter) ** 3)

    # DNA Sequence length
    sequence_length = len(dna_sequence)

    # Longitudinal potential matrix

    global Vz

    Vz = zeros([sequence_length, sequence_length])

    for i in range(sequence_length):
        Vz[i, i] = me * omega[dna_sequence[i], 'zz'] ** 2

    for i in range(sequence_length - 1):
        Vz[i, i + 1] = gamma['zz']
        Vz[i + 1, i] = gamma['zz']

    # Transverse potential matrix

    global Vxy

    Vxy = zeros([2 * sequence_length, 2 * sequence_length])

    for i in range(sequence_length):
        Vxy[i, i] = me * omega[dna_sequence[i], 'xx'] ** 2
        Vxy[sequence_length + i, sequence_length + i] = me * omega[dna_sequence[i], 'yy'] ** 2

    for i in range(sequence_length - 1):
        Vxy[i, i + 1] = gamma['xx']
        Vxy[i + 1, i] = gamma['xx']
        Vxy[sequence_length + i, sequence_length + i + 1] = gamma['yy']
        Vxy[sequence_length + i + 1, sequence_length + i] = gamma['yy']

    for i in range(sequence_length - 1):
        Vxy[i, sequence_length + i + 1] = gamma['xy']
        Vxy[sequence_length + i + 1, i] = gamma['xy']
        Vxy[i + 1, sequence_length + i] = -gamma['xy']
        Vxy[sequence_length + i, i + 1] = -gamma['xy']



# Process eigenvalues and eigenvectors and get relevant quantities

def get_energies_and_amplitudes(eigenvalues, eigenvectors, transverse):
    # Fundamental constants

    ev = 1.602176634e-19 # Electronvolt
    hbar = 6.62607015e-34 / (2 * pi) # Planck's constant
    e_PO = 0.23 # Reference Energy
    me = 9.1093837015e-31 # Electron mass

    # Keep positive eigenvalues only

    negative = []
    for i in range(len(eigenvalues)):
        if eigenvalues[i] < 0:
            negative.append(i)

    values = delete(eigenvalues, negative)
    vectors = delete(eigenvectors, negative, 1)

    # Compute quantities of interest

    energies = sqrt(values / me) * hbar / (2.0 * e_PO * ev)

    if transverse:
        half_comp = int(vectors.shape[0] / 2)
        amplitudes = array([[abs(vectors[comp, vec]) ** 2 + abs(vectors[half_comp + comp, vec]) ** 2 for comp in range(half_comp)] for vec in range(vectors.shape[1])])
    else:
        amplitudes = array([[abs(vectors[comp, vec]) ** 2 for comp in range(vectors.shape[0])] for vec in range(vectors.shape[1])])

    return energies, amplitudes



# Plot energies

def plot_energies():
    # Construct tags and indexes for selected transverse and longitudinal modes

    sel_idx_xy = energies_xy_table.selection
    sel_idx_z = energies_z_table.selection

    if len(sel_idx_xy) + len(sel_idx_z) <= 0:
        return

    sel_idx_xy.sort()
    sel_idx_z.sort()

    tag_xy = [f"T{i}" for i in sel_idx_xy]
    tag_z = [f"L{i}" for i in sel_idx_z]

    tag = tag_xy + tag_z

    # Select frequencies

    sel_energies_xy = energies_xy[sel_idx_xy]
    sel_energies_z = energies_z[sel_idx_z]

    # Construct indexes for colormap

    max_len = max(len(energies_xy) - 1, len(energies_z) - 1)
    norm_idx_xy = [s / max_len for s in sel_idx_xy]
    norm_idx_z = [s / max_len for s in sel_idx_z]
    norm_idx = norm_idx_xy + norm_idx_z

    # Join lists
    freqs = sel_energies_xy.tolist() + sel_energies_z.tolist()

    # Plot options

    fig = Figure(figsize = (plot_width_input.value, plot_height_input.value), layout = 'constrained')
    ax = fig.subplots()

    ax.set_xlabel('Modes')
    ax.set_ylabel('Energy (Epo)')

    color_map = plt.get_cmap(cmaps_picker.value_name)

    ax.bar(tag, freqs, color = color_map(norm_idx))

    # Remove tick labels to enhance visibility

    tick_labels = ax.xaxis.get_majorticklabels()
    if len(tick_labels) >= 20:
        for i, label in enumerate(tick_labels):
            if i % int(len(tick_labels) / 20) > 0:
                label.set_visible(False)

    return fig



# Plot modes

def plot_modes():
    # Construct tags and markers for selected transverse and longitudinal modes

    sel_idx_xy = energies_xy_table.selection
    sel_idx_z = energies_z_table.selection

    if len(sel_idx_xy) + len(sel_idx_z) <= 0:
        return

    sel_idx_xy.sort()
    sel_idx_z.sort()

    tag_xy = [f"T{i}" for i in sel_idx_xy]
    tag_z = [f"L{i}" for i in sel_idx_z]

    tag = tag_xy + tag_z

    marker_xy = ['.-' for i in sel_idx_xy]
    marker_z = ['.:' for i in sel_idx_z]
    marker = marker_xy + marker_z

    # Select modes

    sel_modes_xy = modes_xy[sel_idx_xy]
    sel_modes_z = modes_z[sel_idx_z]

    # Construct indexes for colormap

    max_len = max(len(modes_xy) - 1, len(modes_z) - 1)
    norm_idx_xy = [s / max_len for s in sel_idx_xy]
    norm_idx_z = [s / max_len for s in sel_idx_z]
    norm_seq = norm_idx_xy + norm_idx_z

    # Join lists
    modes = sel_modes_xy.tolist() + sel_modes_z.tolist()

    # Plot options

    fig = Figure(figsize = (plot_width_input.value, plot_height_input.value), layout = 'constrained')
    ax = fig.subplots()

    ax.set_xlabel('DNA Sequence')
    ax.set_ylabel('Amplitude')

    color_map = plt.get_cmap(cmaps_picker.value_name)

    for i, mode in enumerate(modes):
        ax.plot(mode, marker[i], color = color_map(norm_seq[i]), label = tag[i])

    dna_sequence = dna_sequence_input.value.upper()
    nucleotids = [letter for letter in dna_sequence]
    ax.set_xticks(range(len(nucleotids)), nucleotids)

    ax.legend(loc = 'upper left', bbox_to_anchor = (0.05, -0.2), ncol = 10, fontsize = 'x-small')

    return fig



# Energies tables

def set_energies_xy_table_data():
    # Tags

    idx = list(range(len(energies_xy)))
    tag = [f"T{i:03d}" for i in idx]

    # Data

    df = pd.DataFrame({
        'Modes': tag,
        'Energies (Epo)': energies_xy.tolist()
    }, index = idx)

    energies_xy_table.value = df



def set_energies_z_table_data():
    # Tags

    idx = list(range(len(energies_z)))
    tag = [f"L{i:03d}" for i in idx]

    # Data

    df = pd.DataFrame({
        'Modes': tag,
        'Energies (Epo)': energies_z.tolist()
    }, index = idx)

    energies_z_table.value = df



# Set plots

def set_outputs(event):
    fig_energies = plot_energies()
    energies_pane.object = fig_energies

    fig_modes = plot_modes()
    modes_pane.object = fig_modes



pn.bind(set_outputs, plot_button, watch = True)



# Remove plots

def clear_outputs():
    energies_pane.object = None
    modes_pane.object = None



# Main computing function

def compute_vibrational_modes(event):
    dna_sequence = dna_sequence_input.value.upper()

    sequence_length = len(dna_sequence)

    # Checks

    if sequence_length < 2:
        notifications.error('Length of DNA chain smaller than 2!')
        return

    # Build dynamic matrices

    set_dynamic_matrices()

    # Diagonalize dynamic matrices

    global eigenvalues_xy, eigenvectors_xy
    global eigenvalues_z, eigenvectors_z

    eigenvalues_xy, eigenvectors_xy = linalg.eigh(Vxy)
    eigenvalues_z, eigenvectors_z = linalg.eigh(Vz)

    # Get frequencies and modes

    global energies_xy, modes_xy
    global energies_z, modes_z

    energies_xy, modes_xy = get_energies_and_amplitudes(eigenvalues_xy, eigenvectors_xy, True)
    energies_z, modes_z = get_energies_and_amplitudes(eigenvalues_z, eigenvectors_z, False)

    set_energies_xy_table_data()
    set_energies_z_table_data()

    # Clear plots
    clear_outputs()

    # Enable action buttons
    plot_button.disabled = False
    data_download_button.disabled = False



dna_sequence_input.param.watch(compute_vibrational_modes, 'value')
distance_input.param.watch(compute_vibrational_modes, 'value')
twist_angle_input.param.watch(compute_vibrational_modes, 'value')



# Set in-memory zip with input and output data

def get_zip_data():
    # Set data zip filename

    now = datetime.now()
    name = now.strftime("dna_data_%Y-%m-%d_%Hh%Mm%Ss")
    data_download_button.filename = name + '.zip'

    # Prepare parameters data

    parameters_df = pd.DataFrame({
        'DNA Sequence': dna_sequence_input.value.upper(),
        'Distance (Angstrom)': distance_input.value,
        'Twist Angle (rad)': twist_angle_input.value
    }, index = [0])

    # Prepare CSV data

    data_csv = [
        ('parameters.csv', BytesIO(parameters_df.to_csv(index = False, header = True).encode())),
        ('Vxy.csv', BytesIO(pd.DataFrame(Vxy).to_csv(index = False, header = False).encode())),
        ('Vz.csv', BytesIO(pd.DataFrame(Vz).to_csv(index = False, header = False).encode())),
        ('eigenvalues_xy.csv', BytesIO(pd.DataFrame(eigenvalues_xy).to_csv(index = False, header = False).encode())),
        ('eigenvalues_z.csv', BytesIO(pd.DataFrame(eigenvalues_z).to_csv(index = False, header = False).encode())),
        ('eigenvectors_xy.csv', BytesIO(pd.DataFrame(eigenvectors_xy).to_csv(index = False, header = False).encode())),
        ('eigenvectors_z.csv', BytesIO(pd.DataFrame(eigenvectors_z).to_csv(index = False, header = False).encode())),
        ('energies_xy.csv', BytesIO(pd.DataFrame(energies_xy).to_csv(index = False, header = False).encode())),
        ('energies_z.csv', BytesIO(pd.DataFrame(energies_z).to_csv(index = False, header = False).encode())),
        ('amplitudes_xy.csv', BytesIO(pd.DataFrame(modes_xy).to_csv(index = False, header = False).encode())),
        ('amplitudes_z.csv', BytesIO(pd.DataFrame(modes_z).to_csv(index = False, header = False).encode()))
    ]

    # Prepare JSON data

    data_json = [
        ('parameters.json', BytesIO(parameters_df.to_json(orient = 'records').encode())),
        ('Vxy.json', BytesIO(pd.DataFrame(Vxy, index = [f"row {i}" for i in range(Vxy.shape[0])], columns = [f"col {i}" for i in range(Vxy.shape[1])]).to_json(orient = 'index').encode())),
        ('Vz.json', BytesIO(pd.DataFrame(Vz, index = [f"row {i}" for i in range(Vz.shape[0])], columns = [f"col {i}" for i in range(Vz.shape[1])]).to_json(orient = 'index').encode())),
        ('eigenvalues_xy.json', BytesIO(pd.DataFrame([eigenvalues_xy], columns = list(range(len(eigenvalues_xy)))).to_json(orient = 'values').encode())),
        ('eigenvalues_z.json', BytesIO(pd.DataFrame([eigenvalues_z], columns = list(range(len(eigenvalues_z)))).to_json(orient = 'values').encode())),
        ('eigenvectors_xy.json', BytesIO(pd.DataFrame(eigenvectors_xy, index = list(range(eigenvectors_xy.shape[0])), columns = [f"vector {i}" for i in range(eigenvectors_xy.shape[1])]).to_json(orient = 'columns').encode())),
        ('eigenvectors_z.json', BytesIO(pd.DataFrame(eigenvectors_z, index = list(range(eigenvectors_z.shape[0])), columns = [f"vector {i}" for i in range(eigenvectors_z.shape[1])]).to_json(orient = 'columns').encode())),
        ('energies_xy.json', BytesIO(pd.DataFrame([energies_xy], columns = list(range(len(energies_xy)))).to_json(orient = 'values').encode())),
        ('energies_z.json', BytesIO(pd.DataFrame([energies_z], columns = list(range(len(energies_z)))).to_json(orient = 'values').encode())),
        ('amplitudes_xy.json', BytesIO(pd.DataFrame(modes_xy, index = list(range(modes_xy.shape[0])), columns = [f"mode {i}" for i in range(modes_xy.shape[1])]).to_json(orient = 'columns').encode())),
        ('amplitudes_z.json', BytesIO(pd.DataFrame(modes_z, index = list(range(modes_z.shape[0])), columns = [f"mode {i}" for i in range(modes_z.shape[1])]).to_json(orient = 'columns').encode()))
    ]

    # Make zip file to download

    zip_buffer = BytesIO()

    with zipfile.ZipFile(zip_buffer, 'a', zipfile.ZIP_DEFLATED) as archive:
        archive.mkdir(name)
        # CSV
        archive.mkdir(name + '/csv')
        for file_name, file_data in data_csv:
            archive.writestr(name + '/csv/' + file_name, file_data.getvalue())
        # JSON
        archive.mkdir(name + '/json')
        for file_name, file_data in data_json:
            archive.writestr(name + '/json/' + file_name, file_data.getvalue())

    zip_buffer.seek(0)

    return zip_buffer



data_download_button.callback = get_zip_data



# Output layouts

counts_layout = pn.Row(dna_sequence_length_indicator, a_length_indicator, t_length_indicator, c_length_indicator, g_length_indicator, sizing_mode = 'stretch_width')
panes_layout = pn.Column(energies_pane, modes_pane, sizing_mode = 'stretch_both')

# Controls cards

card_margin = 20

computation_parameters_card = pn.Card(distance_input, twist_angle_input, title = 'Parameters', margin = card_margin)
xy_modes_card = pn.Card(energies_xy_table, title = 'Transverse modes', margin = card_margin)
t_modes_card = pn.Card(energies_z_table, title = 'Longitudinal modes', margin = card_margin)
plot_options_card = pn.Card(plot_width_input, plot_height_input, cmaps_picker, title = 'Plot options', margin = card_margin)
actions_card = pn.Card(plot_button, data_download_button, title = 'Actions', margin = card_margin)

controls_layout = pn.FlexBox(xy_modes_card, t_modes_card, computation_parameters_card, plot_options_card, actions_card, flex_direction = 'row', sizing_mode = 'stretch_width')

# Main layout

main_layout = pn.Column(dna_sequence_input, counts_layout,  pn.layout.Divider(), controls_layout, pn.layout.Divider(), panes_layout, sizing_mode = 'stretch_both', scroll = True)

# Template

template = pn.template.BootstrapTemplate(title = 'DNA Electronic Quantum Vibrational Modes', header = links_pane)
template.main.append(main_layout)
template.servable()
