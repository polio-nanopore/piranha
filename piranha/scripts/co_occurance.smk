import os
import collections
from Bio import SeqIO
import yaml

from piranha.analysis.clean_gaps import *
from piranha.analysis.consensus_functions import *
from piranha.utils.log_colours import green,cyan
from piranha.utils.config import *


BARCODE = config[KEY_BARCODE]
SAMPLE = str(config[KEY_SAMPLE])
REFERENCES = config[BARCODE]
