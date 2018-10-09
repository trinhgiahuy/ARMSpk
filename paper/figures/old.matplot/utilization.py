import openpyxl
import urllib
import tempfile
import os

#-- Settings --#

XLSX_LINK = "https://docs.google.com/spreadsheets/d/1un0TIi31LXI9yURmwobPkCOatNXTvVPofXbEu_W-SvA/export?format=xlsx&id=1un0TIi31LXI9yURmwobPkCOatNXTvVPofXbEu_W-SvA"

MAX_COL=100
MAX_ROW=50

TITLE='Utilization HPC Centers'

BOTTOM=-5
COL_SYSTEM=0

def try_float(value):
    try:
        return float(value)
    except ValueError:
        return BOTTOM

AB_TO_LABEL = {
    "geo": "Geoscience/Earthscience",
    "chm": "Chemistry",
    "phy": "Physics",
    "qcd": "Lattice QCD",
    "mat": "Material Science/Engineering",
    "eng": "Engineering (Mechanics, CFD)",
    "mcs": "Math/Computer Science",
    "bio": "Bioscience",
    "oth": "Other"
}

def get_year(row):
    val = search_col(row, "Year")
    if val is None:
        return "XXXX"
    return str(int(val)).replace('20', '\'')

def get_col(row, ablabel):
    label = AB_TO_LABEL[ablabel]
    val = search_col(row, label)
    if val is None:
        return BOTTOM
    return try_float(val)

#--------------#

def temp_dl(dl_link, ext):
    with tempfile.NamedTemporaryFile(delete=False) as temp:
        temp.write(urllib.request.urlopen(dl_link).read())
        loc = temp.name

    new_loc = loc + "." + ext
    os.rename(loc, new_loc)
    return new_loc

XLSX_FILE = temp_dl(XLSX_LINK, "xlsx")
WB = openpyxl.load_workbook(XLSX_FILE, data_only=True)

# RESULT_DIC[NAME OF SYSTEM][COLUMNS]
RESULT_DIC = { # Literal for setting the order of contents
    'ANL': True,
    'NERSC': True,
    'HLRS': True,
    'RRZE': True,
    'CSCS': True,
    'R-CCS K-Computer': True,
    'U. Tokyo Oakforest-PACS': True,
    'NARLabs': True
}

for ws in WB:
    for row in ws.iter_rows(min_row=1, max_col=MAX_COL, max_row=MAX_ROW):
        if row[COL_SYSTEM].value in RESULT_DIC:
            sys = row[COL_SYSTEM].value
            RESULT_DIC[sys] = row

def delete_whitespace(string):
    return ''.join(string.split())

# Without whitespaces
def startswith_wo_ws(a, b):
    return delete_whitespace(a).startswith(delete_whitespace(b))

# Find the position of the columns using the first row of the sheet
def search_col(row, label):
    n = -1
    for ws in WB:
        if ws.title == TITLE:
            for idx, cell in enumerate(next(ws.rows)):
                if not isinstance(cell, type(None)) and startswith_wo_ws(cell.value, label):
                    n = idx
                    break
            if n != -1:
                return row[n].value

    print('NOT FOUND: ' + label + ' of ' + title)
    return None

from draw import *

#figure_setup("Breakdown of System Utilization by Domain", "%", (16, 9), 20)
figure_setup("", "%", (16, 6), 20)

ratio_sys = range(len(RESULT_DIC) * 3)
ratio_label = []
ratio = { "geo": [], "chm": [], "phy": [],
          "qcd": [], "mat": [], "eng": [],
          "mcs":  [], "bio": [], "oth": [] }

for sys, row in RESULT_DIC.items():
    ratio_label.append("")
    if sys == 'R-CCS K-Computer':
        label = 'R-CCS\nK-Computer\n'
    elif sys == 'U. Tokyo Oakforest-PACS':
        label = 'U. Tokyo\nOakforest-\nPACS'
    elif sys == 'NARLabs':
        label = 'NAR-\nLabs'
    else:
        label = sys
    ratio_label.append(label + " (" + get_year(row) + ")")
    ratio_label.append("")

    sum = 0
    for app in list(ratio.keys())[::-1]:
        data = ratio[app]
        data.append(BOTTOM)
        sum+=get_col(row, app)
        data.append(sum)
        data.append(BOTTOM)

plt.xticks(np.arange(len(ratio_sys)), ratio_label)

for app, data in ratio.items():
    plt.bar(ratio_sys, data, 0.8, label=app)

figure_save("img/utilization.eps")
