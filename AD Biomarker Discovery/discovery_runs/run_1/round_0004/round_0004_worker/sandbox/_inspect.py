import sys, os, inspect
print('PYTHONPATH env', os.environ.get('PYTHONPATH'))
print(sys.path[:8])
import shared_analysis
from shared_analysis import sea_ad_lfb
print(shared_analysis.__file__)
print(inspect.signature(sea_ad_lfb.build_cell_table))
