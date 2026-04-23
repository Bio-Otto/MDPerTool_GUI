import io
import logging

from mdpertool.src.logger import Logger, _make_stream_encoding_safe, LOGGER_NAME


def test_make_stream_encoding_safe_handles_none():
    assert _make_stream_encoding_safe(None) is None


def test_make_stream_encoding_safe_keeps_non_reconfigurable_stream():
    stream = io.StringIO()
    assert _make_stream_encoding_safe(stream) is stream


def test_logger_is_idempotent_for_same_file(tmp_path):
    log_file = tmp_path / "run.log"

    logger_first = Logger(str(log_file))
    logger_second = Logger(str(log_file))

    assert logger_first is logger_second
    assert logger_first.name == LOGGER_NAME
    assert len(logger_first.handlers) == 2
    assert all(isinstance(h, logging.Handler) for h in logger_first.handlers)
