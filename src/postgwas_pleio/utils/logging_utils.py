import logging
import sys
from pathlib import Path


# =========================================================
# LOGGER SETUP
# =========================================================
def setup_logger(log_file: Path):

    log_file = Path(log_file)
    log_file.parent.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger("PLEIO_MERGE")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    # -----------------------------------------------------
    # FILE HANDLER → FULL VERBOSE
    # -----------------------------------------------------
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.INFO)     # everything saved
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    # -----------------------------------------------------
    # CONSOLE HANDLER → IMPORTANT ONLY
    # -----------------------------------------------------
    ch = logging.StreamHandler(sys.stderr)
    ch.setLevel(logging.WARNING)  # only warning + error on screen
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    logger.info("=" * 80)
    logger.info("[ACTION] PIPELINE STARTED")
    logger.info("=" * 80)

    return logger


# =========================================================
# SEMANTIC LOGGING HELPERS
# =========================================================
def log_action(logger, msg):
    logger.info(f"[ACTION] {msg}")


def log_detected(logger, msg):
    logger.info(f"[DETECTED] {msg}")


def log_warning(logger, msg):
    logger.warning(f"[WARNING] {msg}")


def log_error(logger, msg):
    logger.error(f"[ERROR] {msg}")