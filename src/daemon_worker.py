#!/usr/bin/env python3
"""
Deamon-like process that checks for jobs and starts analysis on them
"""
import config
from config import logger
import os
from time import sleep, ctime
import threading
from smtplib import SMTP
import sys
import traceback
from subprocess import call


def send_email(to, subject, msg_body):
    """
    Send out an email about the status of the finished job
    """
    msg = f"To: {to}\n"
    msg += f"Subject: {subject}\n\n"
    msg += "Dear DARIO user,\n"
    msg += msg_body

    # Send e-mail
    try:
        server = SMTP()
        server.connect("bierdepot.bioinf.uni-leipzig.de")
        server.sendmail("dario@bioinf.uni-leipzig.de", to, msg)
        server.close()
    except Exception as inst:
        logger.error("Email Address error: " + str(inst))


def workon(job):
    """
    Run on new analysis job
    """
    # Extract job parameters
    cols = job.strip().split("|")
    job_hash = cols[0]
    email = cols[2].strip()
    filename = cols[4]

    # Try to run analyzer
    runpath = config.WORKER_JOBS_PATH + job_hash + "/"
    logger.info(" -- Starting DARIO Analysis in " + runpath + "\n")
    return_code_analysis = call("python " + config.WORKER_SOURCE_PATH + "/analysis.py -i " + runpath, shell=True)
    logger.info(" -- Finished DARIO Analysis in " + runpath + "\n")

    # Update wokon file
    workon_filename = config.WORKER_JOBS_PATH + config.WORKON_JOBS_FILENAME
    with open(workon_filename) as workon_file:
        lines = workon_file.readlines()
    lines.remove(job)
    with open(workon_filename, "w") as workon_file:
        for line in lines:
            workon_file.write(line)

    # Move files to public_html
    return_code = call(["scp", "-r", runpath, f"{config.WEBSERVER_SSH_MACHINE}:{config.WEBSERVER_RESULTS_PATH}"])
    # Workon file on server
    webserver_scp_path = f"{config.WEBSERVER_SSH_MACHINE}:{config.WEBSERVER_JOBS_PATH}"
    call(["scp", workon_filename, webserver_scp_path + config.WORKON_JOBS_FILENAME])

    # Remove files locally
    if return_code == 0:
        return_code = call(["rm", "-r", runpath])
    else:
        logger.error(f"Error during analysis {job_hash}")

    # Write user e-mail
    webpath = f"{config.WEB_URL}sresult/{job_hash}/index.html"
    if email:
        msg = f"The results of your request file {filename} can be found at:\n{webpath}\n"
        send_email(email, "DARIO computation completed", msg)

    # Write internal failure e-mail
    if return_code_analysis != 0 or return_code != 0:
        msg = f"The DARIO job {job_hash} on file {filename} has failed. See:\n{webpath}\n"
        send_email(config.FAILURE_EMAIL, "DARIO job failed", msg)


# The Python main method
if __name__ == "__main__":
    logger.info("The daemon starts")
    while True:
        # list of new jobs
        with open(config.WORKER_JOBS_PATH + config.JOB_QUEUE_FILENAME) as f:
            lines = f.readlines()

        if len(lines) > 0 and threading.activeCount() <= config.WORKER_THREADS:
            # Remove first thread from auto_file
            job = lines[0]
            del lines[0]
            with open(config.WORKER_JOBS_PATH + config.JOB_QUEUE_FILENAME, "w") as f:
                for line in lines:
                    f.write(line)

            # Add to currently-on-work list
            with open(config.WORKER_JOBS_PATH + config.WORKON_JOBS_FILENAME, "a") as workon_file:
                workon_file.write(job)

            # Update Queue files on server
            webserver_scp_path = f"{config.WEBSERVER_SSH_MACHINE}:{config.WEBSERVER_JOBS_PATH}"
            call(
                [
                    "scp",
                    config.WORKER_JOBS_PATH + config.JOB_QUEUE_FILENAME,
                    webserver_scp_path + config.JOB_QUEUE_FILENAME,
                ]
            )
            call(
                [
                    "scp",
                    config.WORKER_JOBS_PATH + config.WORKON_JOBS_FILENAME,
                    webserver_scp_path + config.WORKON_JOBS_FILENAME,
                ]
            )

            # Note: activeCount includes this program
            logger.info(f"Starting Thread: {threading.activeCount()} of {config.WORKER_THREADS}")

            P2 = threading.Thread(target=workon, args=[job])
            P2.start()

            # occasionally reload config to allow change of WORKHORSE_THREADS
            # reload(config)
        sleep(15)

    logger.error("The daemon dies")
