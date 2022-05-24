#!/usr/bin/env python3
import config
from config import logger
from time import sleep
from subprocess import call


def submit_jobs():
    """
    Checks for new jobs that have been received by the webserver, move them
    to the worker machine, and update the jobs status file
    """
    with open(config.WEBSERVER_JOBS_PATH + config.NEW_PROCESS_FILENAME) as f:
        lines = f.readlines()

    # Are there are new processes?
    if len(lines) > 0: 
        # Get job parameters
        job = lines[0].strip()
        cols = job.split("|")
        job_hash = cols[0]

        # Remove first job from new processed file and move job to workhorse queue
        # @todo Beware of race conditions, potentially lock file
        del lines[0]
        with open(config.WEBSERVER_JOBS_PATH + config.NEW_PROCESS_FILENAME, "w") as f:
            for line in lines:
                f.write(line)

        # Move to worker
        rc = 0
        # Copy work directory
        rc = rc | call(
            [
                "scp",
                "-r",
                config.WEBSERVER_JOBS_PATH + job_hash,
                config.WORKER_SSH_MACHINE + ":" + config.WORKER_JOBS_PATH,
            ]
        )

        # Append to job queue...
        rc = rc | call(
            [
                "ssh",
                config.WORKER_SSH_MACHINE,
                "echo '" + job + "'>>" + config.WORKER_JOBS_PATH + config.JOB_QUEUE_FILENAME,
            ]
        )

        # ... and update local mirror of job list
        rc = rc | call(
            [
                "scp",
                config.WORKER_SSH_MACHINE + ":" + config.WORKER_JOBS_PATH + config.JOB_QUEUE_FILENAME,
                config.WEBSERVER_JOBS_PATH,
            ]
        )

        # @todo give back some error to the user if return_code != 0
        logger.info(f"Job {job_hash} has been added")


if __name__ == "__main__":
    logger.info("The daemon starts")
    try:
        while True:
            submit_jobs()
            sleep(15)
    except:
        logger.error("The daemon dies")
