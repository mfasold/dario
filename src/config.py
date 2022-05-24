"""
Sets the main config variables used in DARIO    
"""

# Set up logger
import logging

logging.basicConfig(format="%(levelname)s:%(asctime)s %(message)s", level=logging.INFO)
logger = logging.getLogger("DARIO")


WEBSERVER_SSH_MACHINE = "maitai2"
WEBSERVER_SOURCE_PATH = "/u/dario/dario-app/src/"
WEBSERVER_TEMPLATES_PATH = "/u/dario/dario-app/src/templates/"
WEBSERVER_JOBS_PATH = "/u/dario/computations/wrk/"
WEBSERVER_RESULTS_PATH = "/u/dario/public_html/result/"
WEBSERVER_EXAMPLE_DATA_PATH = "/u/dario/public_html/example/"

WORKER_SSH_MACHINE = "k74"
WORKER_TEMPLATES_PATH = "/scratch/dario/src/templates/"
WORKER_JOBS_PATH = "/scratch/dario/computations/"
WORKER_SOURCE_PATH = "/scratch/dario/src/"
WORKER_THREADS = 3

WEB_URL = "http://dario.bioinf.uni-leipzig.de/"

NEW_PROCESS_FILENAME = "new_process.list"
JOB_QUEUE_FILENAME = "auto.list"
ALL_JOBS_FILENAME = "all_jobs.list"
WORKON_JOBS_FILENAME = "workon.list"
PARAMS_FILENAME = "job_params.txt"
STDERR_FILENAME = "stderror.log"
RUNLOG_FILENAME = "run2.log"
ERROR_FILENAME = "error.log"

# Paramters used by cronjob on WORKHORSE compter
WORKHORSE_HOME = "/scratch/dario/"
WORKHORSE_WRKPATH = WORKHORSE_HOME + "computations/"
WORKHORSE_SRCPATH = WORKHORSE_HOME + "src/"
WORKHORSE_THREADS = 3

WEBPATH = "http://dario.bioinf.uni-leipzig.de/"

# Supported annotations
ANNOTATION_DIR_WH = "/scratch/dario/data/annotations/"
SPECIES = {
    "Human (hg18)": {"dir": ANNOTATION_DIR_WH + "hg18/", "id": "hg18", "test_data": "GSM450599.hg18.bed.gz"},
    "Human (hg19)": {"dir": ANNOTATION_DIR_WH + "hg19/", "id": "hg19"},
    "Worm (ce6)": {"dir": ANNOTATION_DIR_WH + "ce6/", "id": "ce6", "test_data": "ce6.GSE17153.bed.gz"},
    "Fruit Fly (dm3)": {"dir": ANNOTATION_DIR_WH + "dm3/", "id": "dm3", "test_data": "dm3.GSE17153.bed.gz"},
    "Rhesus Monkey (rheMac2)": {
        "dir": ANNOTATION_DIR_WH + "rheMac2/",
        "id": "rheMac2",
        "test_data": "rheMac2.GSM450611.bed.gz",
    },
    "Mouse (mm9)": {"dir": ANNOTATION_DIR_WH + "mm9/", "id": "mm9", "test_data": "mm9.GSM314552.bed.gz"},
    "Mouse (mm10)": {"dir": ANNOTATION_DIR_WH + "mm10/", "id": "mm10"},
    "Zebrafish (danRer6)": {"dir": ANNOTATION_DIR_WH + "danRer6/", "id": "danRer6"},
    "Spotted Garr (lepOcu1)": {"dir": ANNOTATION_DIR_WH + "lepOcu1/", "id": "lepOcu1"},
}

RETURN_CODE_OK = 0
RETURN_CODE_ERROR = 1

FAILURE_EMAIL = "dario@bioinf.uni-leipzig.de"
