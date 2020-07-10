import pytest
import common
from testfixtures import log_capture


@log_capture()
@pytest.mark.parametrize(
    "input_file,output_file",
    [('1FAS.pdb', '1FAS_pdb.pqr'), ('1A1P.pdb', '1A1P_assign-only_whitespace_ff=AMBER.pqr')]
)
def test_log_output_in_pqr_location(capture, input_file, output_file):
    """Test that a log file is generated when calling the __init__.py"""
    args = "--log-level=INFO --ff=AMBER"
    tmp_path = common.DATA_DIR
    input_path = common.DATA_DIR / input_file
    output_pqr = output_file
    common.run_pdb2pqr(
        args=args,
        input_pdb=input_path,
        output_pqr=output_pqr,
        tmp_path=tmp_path
    )

    # This match tests that the log file is output to the same location as the output pqr
    capture.check_present(
        ('pdb2pqr.input_output', 'INFO', f'Logs stored: tests/data/{output_file.split(".")[0]}.log')
    )