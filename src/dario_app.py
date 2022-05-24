#!/usr/bin/env python3
import os
import sys
from subprocess import call
from mako.template import Template
from mako.lookup import TemplateLookup
from hashlib import md5
import traceback
from time import ctime, localtime, strftime
from itertools import groupby


from flask import Flask, request, stream_with_context, send_from_directory

app = Flask(__name__)

sys.path.append(os.path.dirname(app.instance_path))  # to import config.py file from the same path
import config as config

app.config["MAX_CONTENT_LENGTH"] = 70 * 1024 * 1024


def render_template(template_name, **kwargs):
    """
    Custom rendering function that uses the Mako template engine instead of Flask-default Jinja2
    """
    mako_lookup = TemplateLookup(directories=[config.WEBSERVER_TEMPLATES_PATH])
    mytemplate = Template(filename=config.WEBSERVER_TEMPLATES_PATH + template_name + ".mako.html", lookup=mako_lookup)
    return mytemplate.render(WEBPATH=config.WEB_URL, phase=3, **kwargs)


@app.route("/")
@app.route("/index.py")  # legacy link compatibility
def index():
    return render_template("index")


@app.route("/help")
@app.route("/privacy")
@app.route("/test_data")
@app.route("/new_job")
def render_static():
    """
    Returns a rendered mako site based on the files in the the templates folder
    """
    return render_template(request.path[1:])


@app.route("/result/<path:path>")
def serve_result(path):
    return send_from_directory(config.WEBSERVER_RESULTS_PATH, path)


@app.route("/example/<path:path>")
def serve_example(path):
    return send_from_directory(config.WEBSERVER_EXAMPLE_DATA_PATH, path)


@app.route("/news")
def render_news():
    """
    Render the news page
    """

    def fasta_iter(fasta_name):
        """
        Given a fasta file, yield tuples of header, sequence (from https://www.biostars.org/p/710/)
        """
        with open(fasta_name) as f:
            faiter = (x[1] for x in groupby(f, lambda line: line[0] == ">"))
            for header in faiter:
                header = header.__next__()[1:].strip()  # drop the ">"
                seq = "".join(s for s in faiter.__next__())
                yield header, seq

    blog_entries = fasta_iter(config.WEBSERVER_TEMPLATES_PATH + "news.fasta")
    return render_template("blog", entries=blog_entries)


@app.route("/upload", methods=["POST"])
def upload_file():
    """
    Handles the upload of data and adding analysis to queue. Implemented
    a "poor-mans" progress, which sends JS responses to the clients in a
    streaming fashion.
    @improvement: Do a proper file upload using a JS library as Dropzone.
    """

    def streamed_response():
        yield render_template("upload")
        try:
            # Get variables & checking
            form = request.form
            use_test_data = False
            if form.get("use_test_data"):
                use_test_data = True

            email = form.get("email")
            species_code = form.get("code")
            coverage_file = request.files["coverage_file"]

            if not use_test_data:
                filename = coverage_file.filename
            else:
                filename = config.SPECIES[species_code]["test_data"]

            # Create hash
            hash = md5((email + filename + ctime()).encode("utf-8")).hexdigest()

            # Begin with upload
            yield '<script type="text/javascript">document.getElementById("Begin").style.visibility = "visible";</script>'

            dir = "%s%s/" % (config.WEBSERVER_JOBS_PATH, hash)
            os.mkdir(dir)

            if not use_test_data:
                coverage_file.save(dir + "mapping_loci.upload")
            else:
                call(["cp", config.WEBSERVER_EXAMPLE_DATA_PATH + filename, dir + "mapping_loci.upload"])

            # Upload user_annotation file
            user_annotation = "NONE"
            if not use_test_data and "user_annotation" in request.files:
                user_annotation_file = request.files["user_annotation"]
                user_annotation_file.save(dir + "mapping_loci.upload")

            # Store run parameters
            jtime = strftime("%Y-%m-%d %H:%M:%S", localtime())
            params = {
                "hash": hash,
                "email": email,
                "code": species_code,
                "filename": filename,
                "job_received_at": jtime,
                "total_upload_size": request.content_length,
                "user_annotation": user_annotation,
                "use_test_data": use_test_data,
            }
            with open(dir + config.PARAMS_FILENAME, "w") as f:
                for k, v in params.items():
                    f.write(k + "\t" + str(v) + "\n")

            # Add to queue
            with open(config.WEBSERVER_JOBS_PATH + config.NEW_PROCESS_FILENAME, "a") as f:
                f.write("%s|NONE|%s|%s|%s\n" % (hash, email, species_code, filename))

            # Add to debug file which stores all started jobs
            remote_ip = request.remote_addr[:-3] + "0"
            with open(config.WEBSERVER_JOBS_PATH + config.ALL_JOBS_FILENAME, "a") as f:
                f.write("%s|%s|%s|%s|%s|%s\n" % (jtime, hash, species_code, filename, email, remote_ip))

            # Now finish upload and tell user
            yield '<script type="text/javascript">document.getElementById("Finish").style.visibility = "visible";</script>'
            yield '<script type="text/javascript">window.setTimeout("location.replace(\\"wait/%s\\")", 3000);</script>' % (
                hash
            )
        except:
            yield '<script type="text/javascript">alert("Exception!");</script>'

    return app.response_class(stream_with_context(streamed_response()))


@app.errorhandler(413)
def upload_too_large(e):
    return "You file exceeded maximum size of 60MB. Please try to compress your file, using e.g. gzip.", 413


@app.route("/wait/<job_hash>")
def show_wait_page(job_hash):
    """
    Shows the user his position in the job queue
    """
    found_msg = ""
    with open(config.WEBSERVER_JOBS_PATH + config.JOB_QUEUE_FILENAME) as f:
        for i, line in enumerate(f):
            if line.split("|", 1)[0] == job_hash:
                found_msg = f"Your job is on position {i+1} in the queue."
                break
    if found_msg == "":
        with open(config.WEBSERVER_JOBS_PATH + config.WORKON_JOBS_FILENAME) as f:
            for line in f:
                if line.split("|", 1)[0] == job_hash:
                    found_msg = "Your job is on work.<br>\n"
                    break
    if found_msg == "" and os.path.exists(config.WEBSERVER_RESULTS_PATH + job_hash + "/index.html"):
        found_msg = "Your job is finished. You will forwarded to the results."
        found_msg += (
            f'<script type="text/javascript">window.setTimeout("location.replace(\\"/result/%s/index.html\\")", 5000);</script>'
            % (job_hash)
        )
    if found_msg == "":
        found_msg = "Adding your job to the queue."
    return render_template("job_wait", found=found_msg, pagelink=f"wait/{job_hash}")


@app.route("/finished_jobs/<secret_code>")
def show_finished_jobs(secret_code):
    """
    Provides a page with finished jobs
    """
    if 1:
        with open(config.WEBSERVER_JOBS_PATH + config.ALL_JOBS_FILENAME) as f:
            lines = [x.split("|") for x in f.readlines()]
        for job in lines:
            job.append(os.path.exists(config.WEBSERVER_RESULTS_PATH + job[1] + "/ncRNA.expression.bed"))
            job.append(os.path.exists(config.WEBSERVER_RESULTS_PATH + job[1] + "/predictions.expression.bed"))
        return render_template("finished_jobs", job_infos=lines)
    else:
        return "Access not permitted"


application = app  # @todo: remove. WSGI script reloading for testing
