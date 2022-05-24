#!/usr/bin/python
import os
import sys
from subprocess import check_call, call, PIPE, STDOUT
from time import localtime, strftime
import config
from config import logger

SRCDIR = config.WORKER_SOURCE_PATH

# for HTML generation
import csv
from mako.template import Template
from mako.lookup import TemplateLookup


def render_template(template_name, **kwargs):
    """
    Custom rendering function that uses the Mako template engine (instead of Flask-default Jinja2)
    """
    mako_lookup = TemplateLookup(directories=[config.WORKER_TEMPLATES_PATH])
    mytemplate = Template(filename=config.WORKER_TEMPLATES_PATH + template_name + ".mako.html", lookup=mako_lookup)
    return mytemplate.render(WEBPATH=config.WEB_URL, phase=2, **kwargs)


def create_expression_table_HTML(
    rundir,
    table_file="",
    template_file="result_expression_table",
    html_prefix="ncRNA_table_",
    sort_key=lambda n: n[4],
    reverse=True,
):
    """
    Create a dict containing, for each ncRNA type, lists of file rows
    """
    ncrna_expression = {}
    ncrna_expression_file = csv.reader(open(rundir + table_file, "rU"), delimiter="\t")
    for row in ncrna_expression_file:
        ncrna_type = row[6]
        if ncrna_type in ncrna_expression:
            ncrna_expression[ncrna_type].append(row)
        else:
            ncrna_expression[ncrna_type] = [row]

    # Sort rows according to expression
    for ncrna_type in list(ncrna_expression.keys()):
        ncrna_expression[ncrna_type] = sorted(ncrna_expression[ncrna_type], key=sort_key, reverse=reverse)

    # Render HTML file
    for ncrna_type in list(ncrna_expression.keys()):
        with open(rundir + html_prefix + ncrna_type + ".html", "w") as f:
            f.write(
                render_template(
                    template_file, ncrna_expression=ncrna_expression[ncrna_type], ncrna_type=ncrna_type, rundir=rundir
                )
            )

    # Return ncRNA types and number of items in each class
    return dict((k, len(ncrna_expression[k])) for k in list(ncrna_expression.keys()))


def swap_ncrna_to_index(seq, key, new_index):
    """Swaps list elements that new position of element e with e[0] == key is new_index"""
    el_index = next((i for i in range(len(seq)) if seq[i][0] == key), None)
    if el_index >= 0 and new_index < len(seq):  # only if found and possible
        seq[new_index], seq[el_index] = seq[el_index], seq[new_index]  # swap
        return True
    return False


def sort_ncrna_records(seq, order_list):
    """Orders list of records r such that r[0] is ordered according to order_list"""
    current_top_index = 0
    for ncrna in order_list:
        if swap_ncrna_to_index(seq, ncrna, current_top_index):
            current_top_index += 1


def get_classifier_statistics(statfile="my.modelstat"):
    """Parse classifier information out of my.modelstat"""
    lines = [line.strip() for line in open(statfile)]

    # Remove anything before cross validation
    lines = lines[lines.index("=== Stratified cross-validation ===") :]

    # Now parse textfile 'dirty'
    matrix_idx = lines.index("=== Confusion Matrix ===") + 3
    ncrna_count = 4
    confusion_matrix = []
    ncrna_classes = []
    for k in range(0, ncrna_count):
        line_split = lines[matrix_idx + k].split()
        confusion_matrix.append(line_split[0:ncrna_count])
        ncrna_classes.append(line_split[-1])
    return {"ncrna_count": ncrna_count, "ncrna_classes": ncrna_classes, "confusion_matrix": confusion_matrix}


def create_main_HTML(rundir, job_infos):
    """Creates Result HTML"""
    # Get data needed for the template
    analysis_infos = {}
    info_file = csv.reader(open(rundir + "upload.info", "rU"), delimiter=":")
    for row in info_file:
        if len(row) > 1:
            analysis_infos[row[0]] = row[2].strip()

    analysis_infos["results_dir"] = config.WEB_URL + "result/" + job_infos["hash"]

    # Create list containg class, #genes, reads and expression for each ncRNA class
    reads_infos = {}
    reads_file = csv.reader(open(rundir + "reads.info", "rU"), delimiter=":")
    for row in reads_file:
        if len(row) > 1:
            reads_infos[row[0]] = [row[3].strip(), row[4].strip()]

    # Fetch data for ncRNA expression
    ncrna_expression = create_expression_table_HTML(
        rundir, "ncRNA.expression.bed", sort_key=lambda n: n[0] + "a" + n[1], reverse=False
    )
    analysis_infos["ncrna_expression"] = [
        [k, genes, reads_infos[k][0], reads_infos[k][1]] for k, genes in list(ncrna_expression.items())
    ]
    sort_ncrna_records(analysis_infos["ncrna_expression"], ["miRNA", "snoRNA_CD", "snoRNA_HACA", "tRNA"])

    # Fetch data for predictions
    if job_infos["prediction_succesful"] == "1":
        predictions = create_expression_table_HTML(
            rundir, "predictions.expression.bed", "result_predictions_table", "predictions_table_"
        )
        analysis_infos["predictions"] = [[k, genes, 0, 0] for k, genes in list(predictions.items())]
        sort_ncrna_records(analysis_infos["predictions"], ["miRNA", "snoRNA_CD", "snoRNA_HACA", "tRNA"])

        # Get classifier statistics, if possible
        try:
            analysis_infos["classifier"] = get_classifier_statistics(rundir + "my.modelstat")
        except:
            pass

    # Fetch data for user annotation - similar to ncRNA expression above
    if job_infos["user_annotation"] != "NONE" and job_infos["user_annotation_succesful"] == "1":
        user_annotation = create_expression_table_HTML(
            rundir,
            "user_annotation.expression.bed",
            html_prefix="user_annotation_table_",
            sort_key=lambda n: n[0] + "a" + n[1],
            reverse=False,
        )
        analysis_infos["user_annotation"] = [[k, genes] for k, genes in list(user_annotation.items())]

    # Render and write template
    with open(rundir + "index.html", "w") as f:
        f.write(render_template("result_index", analysis=analysis_infos, job_infos=job_infos, rundir=rundir))


def create_error_HTML(rundir, job_infos, error_msg="", error_pre=""):
    """"Create result package with an error message"""
    with open(rundir + "index.html", "w") as f:
        f.write(
            render_template("job_error", job_infos=job_infos, rundir=rundir, error_msg=error_msg, error_pre=error_pre)
        )


def try_to_extract_uploaded_file(filename, target_file_basename, original_filename, rundir, stderr):
    """Tries to decompress uploaded file"""
    # Save a zipped copy of origignal filename (may get removed later)
    check_call("gzip -c " + rundir + filename + " >" + rundir + filename + ".gz", shell=True, stderr=stderr)

    # Extract
    supported_archives = (".zip", ".gz", ".tar.gz")
    if original_filename.lower().endswith(supported_archives):  # if archive
        logger.info("Trying to extract " + original_filename)
        check_call(
            ["/bin/bash", SRCDIR + "analysis/extract.sh", rundir, filename, original_filename, target_file_basename],
            stderr=stderr,
        )
    else:  # not an archive -> just copy uploaded file to target name + original file suffix
        logger.info("Copying file" + original_filename)
        original_filename_suffix = os.path.splitext(original_filename)[1].lower()  # get file extension
        check_call(["cp", rundir + filename, rundir + target_file_basename + original_filename_suffix])

    # remove (potentially uncompressed) original to save space (unchecked, might be removed by gunzip earlier)
    call(["rm", rundir + filename], stderr=stderr)


def analyze(rundir):
    """Main driver script"""
    if not rundir.endswith("/"):  # for robustness, add "/"
        rundir = rundir + "/"

    # Setup error log file
    stderr = open(rundir + config.STDERR_FILENAME, "aw")
    runlog_file = open(rundir + config.RUNLOG_FILENAME, "aw")

    job_infos = {}
    info_file = csv.reader(open(rundir + config.PARAMS_FILENAME, "rU"), delimiter="\t")
    for row in info_file:
        if len(row) > 1:
            job_infos[row[0]] = row[1].strip()

    job_infos["job_start_time"] = strftime("%Y-%m-%d %H:%M:%S", localtime())

    # ------------------------------------------------------------------------------------------
    # STEP 1: Make sure BED or BAM file is present. Extract compressed archive, if necessary
    # ------------------------------------------------------------------------------------------
    mapping_loci_basename = "mapping_loci"
    try:
        try_to_extract_uploaded_file(
            "mapping_loci.upload", mapping_loci_basename, job_infos["filename"], rundir, stderr
        )

        if os.path.exists(rundir + mapping_loci_basename + ".bed"):
            job_infos["mapping_loci_filetype"] = "bed"
        elif os.path.exists(rundir + mapping_loci_basename + ".bam"):
            job_infos["mapping_loci_filetype"] = "bam"
        else:
            raise IOError("No bam/bed file found")

        logger.info(job_infos["mapping_loci_filetype"] + " file found.")
    except:
        stderr.close()
        stderr_text = "\n".join([line.strip() for line in open(rundir + config.STDERR_FILENAME)])
        create_error_HTML(
            rundir,
            job_infos,
            "The uploaded file or archive did not contain or contain more than the required BED or BAM file. Ideally, the filename suffix should indicate its type and therefore end with either .bed or .bam.",
            stderr_text,
        )
        return config.RETURN_CODE_ERROR

    # ------------------------------------------------------------------------------------------
    # STEP 2: Extract user_annotation
    # ------------------------------------------------------------------------------------------
    if job_infos["user_annotation"] != "NONE":
        try:
            try_to_extract_uploaded_file(
                "user_annotation.upload", "user_annotation", job_infos["user_annotation"], rundir, stderr
            )
            if not os.path.exists(rundir + "user_annotation.bed"):
                raise IOError("No user_annotation bam/bed file found")
        except:
            create_error_HTML(rundir, job_infos, "User annotation could not be extracted")
            return config.RETURN_CODE_ERROR

    # ------------------------------------------------------------------------------------------
    # STEP 3: Convert BAM to BED if necessary
    # ------------------------------------------------------------------------------------------
    if job_infos["mapping_loci_filetype"] == "bam":
        try:
            check_call(
                [
                    "perl",
                    SRCDIR + "analysis/map2bed.pl",
                    "-i",
                    rundir + mapping_loci_basename + ".bam",
                    "-f",
                    "3",
                    "-o",
                    rundir + mapping_loci_basename + ".bed",
                    "-z",
                ],
                stdout=runlog_file,
                stderr=stderr,
            )

            # Remove uneeded BAM files
            check_call(["rm", rundir + mapping_loci_basename + ".bam"], stderr=stderr)

            if not os.path.exists(rundir + mapping_loci_basename + ".bed"):
                raise IOError("BED file not present")
        except:
            stderr.close()
            stderr_text = "\n".join([line.strip() for line in open(rundir + config.STDERR_FILENAME)])
            create_error_HTML(rundir, job_infos, "BAM file could not be converted properly.", stderr_text)
            return config.RETURN_CODE_ERROR
        logger.info("Converted BAM 2 BED format")

    # ------------------------------------------------------------------------------------------
    # STEP 4: Check integrity of BED FILE
    # ------------------------------------------------------------------------------------------
    mapping_loci_filename = rundir + mapping_loci_basename + ".bed"
    try:
        # ./checkBed.pl -i reads.bed
        check_call([SRCDIR + "analysis/checkBed.pl", "-i", mapping_loci_filename], stdout=runlog_file, stderr=stderr)
        # ./readsToTags.pl -i reads.bed -s upload.info -l length.out -m multipleMappings.out -o tags.bed
        check_call(
            [
                SRCDIR + "analysis/readsToTags.pl",
                "-i",
                mapping_loci_filename,
                "-s",
                rundir + "upload.info",
                "-l",
                rundir + "length.out",
                "-m",
                rundir + "multipleMappings.out",
                "-o",
                rundir + "tags.bed",
            ],
            stdout=runlog_file,
            stderr=stderr,
        )
    except:
        stderr.close()
        runlog_file.close()
        stderr_text = "<br\>".join([line.strip() for line in open(rundir + config.STDERR_FILENAME)])
        runlog_text = "\n".join([line.strip() for line in open(rundir + config.RUNLOG_FILENAME)])

        if stderr_text.find("memory") != -1:
            create_error_HTML(
                rundir,
                job_infos,
                "We currently do not have sufficent memory to process your file in acceptable time. Sorry, we are working on getting better machines!",
            )
        else:
            create_error_HTML(rundir, job_infos, "Your mapping file has invalid file format.", runlog_text)
        return config.RETURN_CODE_ERROR

    # ------------------------------------------------------------------------------------------
    # REWRITE DONE UNTIL HERE. @todo: Rename files accordingly in the remaining.
    # ------------------------------------------------------------------------------------------

    call(["mv", rundir + "tags.bed", rundir + "upload.bed"])  # rename to upload.bed

    # ------------------------------------------------------------------------------------------
    # STEP 5: Run quantification
    # ------------------------------------------------------------------------------------------
    # annotated ncRNAs
    species = config.SPECIES[job_infos["code"]]
    try:
        check_call(
            [
                "Rscript",
                SRCDIR + "analysis/overlap.R",
                species["dir"] + "ncRNAs.bed",
                species["dir"] + "exons.bed",
                species["dir"] + "introns.bed",
                rundir + "upload.bed",
                rundir + "ncRNAs.reads",
                rundir + "unknown.reads",
                rundir + "reads.info",
                rundir,
            ],
            stderr=stderr,
        )

        check_call(
            [
                "perl",
                SRCDIR + "analysis/writeWig.pl",
                "-i",
                rundir + "ncRNAs.reads",
                "-o",
                rundir,
                "-p",
                "ncRNAs.pos.wig",
                "-n",
                "ncRNAs.neg.wig",
            ],
            stderr=stderr,
        )
        check_call(
            [
                "perl",
                SRCDIR + "analysis/getExpression.pl",
                "-b",
                rundir + "ncRNAs.reads",
                "-o",
                rundir,
                "-a",
                species["dir"] + "ncRNAs.bed",
                "-f",
                rundir + "ncRNA.expression.bed",
                "-p",
                "ncRNAs.pos.wig",
                "-n",
                "ncRNAs.neg.wig",
                "-s",
                species["id"],
            ],
            stderr=stderr,
        )

        ## do some renaming for arabidopsis
        # if species['id'].startswith("ath"):
        #    call("mv {0}ncRNAs.pos.wig {0}ncRNAs.pos.wig.save".format(rundir), shell=True, stderr = stderr)
        #    call("mv {0}ncRNAs.neg.wig {0}ncRNAs.neg.wig.save".format(rundir), shell=True, stderr = stderr)
        #    call("sed 's/=chr/=Chr/' {0}ncRNAs.pos.wig.save | sed '/chloroplast/q' | head -n -1 > {0}ncRNAs.pos.wig".format(rundir), shell=True, stderr = stderr)
        #    call("sed 's/=chr/=Chr/' {0}ncRNAs.neg.wig.save | sed '/chloroplast/q' | head -n -1 > {0}ncRNAs.neg.wig".format(rundir), shell=True, stderr = stderr)
    except:
        create_error_HTML(rundir, job_infos, "Your reads could not be overlapped with ncRNA annotations.")
        return config.RETURN_CODE_ERROR

    # ------------------------------------------------------------------------------------------
    # STEP 6: Run prediction of new candidate ncRNAs
    # ------------------------------------------------------------------------------------------
    try:
        check_call(
            ["/bin/bash", SRCDIR + "analysis/runblockbuster.sh", SRCDIR, rundir, species["dir"] + "ncRNAs.bed"],
            stderr=stderr,
        )
        check_call(
            [
                "perl",
                SRCDIR + "analysis/trainClassifier.pl",
                "-trainSet",
                rundir + "ncRNAs.clusters.flagged",
                "-model",
                rundir + "my.model",
                "-statistic",
                rundir + "my.modelstat",
                "-o",
                rundir,
            ],
            stdout=stderr,
            stderr=STDOUT,
        )
        if not os.path.exists(rundir + "my.model"):
            check_call(["cp", species["dir"] + "my.model", rundir + "my.model"])
        check_call(
            [
                "perl",
                SRCDIR + "analysis/runClassifier.pl",
                "-testSet",
                rundir + "unknown.clusters",
                "-model",
                rundir + "my.model",
                "-predictions",
                rundir + "predictions.bed",
                "-o",
                rundir,
            ],
            stdout=stderr,
            stderr=STDOUT,
        )

        # Quantify predictions
        check_call(
            [
                "Rscript",
                SRCDIR + "analysis/overlapPredictions.R",
                rundir + "predictions.bed",
                rundir + "unknown.reads",
                rundir + "predictions.reads",
                rundir,
            ],
            stderr=stderr,
        )
        check_call(
            [
                "perl",
                SRCDIR + "analysis/writeWig.pl",
                "-i",
                rundir + "predictions.reads",
                "-o",
                rundir,
                "-p",
                "predictions.pos.wig",
                "-n",
                "predictions.neg.wig",
            ],
            stderr=stderr,
        )
        check_call(
            [
                "perl",
                SRCDIR + "analysis/getExpression.pl",
                "-b",
                rundir + "predictions.reads",
                "-o",
                rundir,
                "-a",
                rundir + "predictions.bed",
                "-f",
                rundir + "predictions.expression.bed",
                "-p",
                "predictions.pos.wig",
                "-n",
                "predictions.neg.wig",
                "-r",
                species["dir"] + "RNAz.bed",
                "-s",
                species["id"],
            ],
            stderr=stderr,
        )
        job_infos["prediction_succesful"] = "1"
    except:
        job_infos["prediction_succesful"] = "0"

    # ------------------------------------------------------------------------------------------
    # STEP 7: Quantify user annotation
    # ------------------------------------------------------------------------------------------
    job_infos["user_annotation_succesful"] = "0"  # must be set always # TODO: no! alter create main HTML
    if job_infos["user_annotation"] != "NONE":
        try:
            check_call(
                [
                    "Rscript",
                    SRCDIR + "analysis/overlapPredictions.R",
                    rundir + "user_annotation.bed",
                    rundir + "upload.bed",
                    rundir + "user_annotation.reads",
                ],
                stderr=stderr,
            )
            check_call(
                [
                    "perl",
                    SRCDIR + "analysis/writeWig.pl",
                    "-i",
                    rundir + "user_annotation.reads",
                    "-o",
                    rundir,
                    "-p",
                    "user_annotation.pos.wig",
                    "-n",
                    "user_annotation.neg.wig",
                ],
                stderr=stderr,
            )
            check_call(
                [
                    "perl",
                    SRCDIR + "analysis/getExpression.pl",
                    "-b",
                    rundir + "user_annotation.reads",
                    "-o",
                    rundir,
                    "-a",
                    rundir + "user_annotation.bed",
                    "-f",
                    rundir + "user_annotation.expression.bed",
                    "-p",
                    "user_annotation.pos.wig",
                    "-n",
                    "user_annotation.neg.wig",
                    "-r",
                    species["dir"] + "RNAz.bed",
                    "-s",
                    species["id"],
                ],
                stderr=stderr,
            )
            job_infos["user_annotation_succesful"] = "1"
        except:
            job_infos["user_annotation_succesful"] = "0"  # do nothing
            logger.exception("Exception during analysis")
            # create_error_HTML(rundir, job_infos, "The user annotation could not be processed." )
            # return config.RETURN_CODE_ERROR

    # ------------------------------------------------------------------------------------------
    # STEP 7: Create pictures,result HTML and cleanup
    # ------------------------------------------------------------------------------------------
    job_infos["job_finish_time"] = strftime("%Y-%m-%d %H:%M:%S", localtime())

    try:
        # Create nice figures
        check_call(["Rscript", SRCDIR + "analysis/createQualityControlFigures.R", rundir], stderr=stderr)

        create_main_HTML(rundir, job_infos)
        job_infos["job_completed_succesfully"] = 1
    except:
        create_error_HTML(rundir, job_infos, "Quality control and analysis statistics could not be created.")
        logger.exception("Exception during analysis")
        return config.RETURN_CODE_ERROR

    # Zip all files into one
    call(
        "zip -j {0}dario_results {0}ncRNA.expression.bed {0}predictions.expression.bed {0}user_annotation.expression.bed {0}*.big.png".format(
            rundir
        ),
        shell=True,
        stderr=stderr,
    )

    # Remove unnecessary files
    files_to_delete = [
        mapping_loci_basename + ".bed",
        "upload.bed",
        "unknown.reads",
        "ncRNAs.reads",
        "unknown.clusters",
        "predictions.reads",
        "ncRNAs.clusters.flagged",
        "ncRNAs.clusters",
    ]
    call("rm " + " ".join([rundir + f for f in files_to_delete]), shell=True, stderr=stderr)

    # Write updated job infos file
    with open(rundir + config.PARAMS_FILENAME, "w") as f:
        for k, v in list(job_infos.items()):
            f.write(k + "\t" + str(v) + "\n")

    stderr.close()
    runlog_file.close()

    return config.RETURN_CODE_OK


# Define main entrypoint, to allow usage as standalone tool
if __name__ == "__main__":
    from optparse import OptionParser

    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option("-i", "--input", dest="dir", action="store", type="str", help="The Input Directory")
    (options, args) = parser.parse_args()
    analyze(options.dir)
