"""Input compatibility checks, not a complete FORTRAN namelist-emulator port."""

from pathlib import Path

from pydatcom.io.namelist_parser import NamelistParser


def _flight(text):
    return NamelistParser().parse(text)[0]['namelists']['FLTCON']


def test_unindexed_continuation_preserves_scalars_and_stops_at_assignment():
    result = _flight('$FLTCON NMACH=3., MACH=.6,.9,\n1.4, LOOP=1., HYPERS=.FALSE.$')
    assert result == {'NMACH': 3., 'MACH': [.6, .9, 1.4], 'LOOP': 1., 'HYPERS': False}


def test_indexed_continuation_and_reassignment():
    result = _flight('$FLTCON MACH=.6,.9,1.4,2.5, MACH(2)=.8,1.2, MACH=0.5$')
    assert result['MACH'] == [.5, .8, 1.2, 2.5]


def test_indexed_assignment_promotes_preceding_single_value():
    result = _flight('$FLTCON MACH=.6, MACH(3)=1.4,2.5$')
    assert result['MACH'] == [.6, None, 1.4, 2.5]


def test_repeated_blocks_merge_indexed_and_unindexed_updates():
    result = _flight('$FLTCON MACH=.6,.9,1.4,2.5$\n'
                     '$FLTCON MACH(2)=.8,1.2$\n$FLTCON MACH=.5$')
    assert result['MACH'] == [.5, .8, 1.2, 2.5]
    result = _flight('$FLTCON MACH=.6$\n$FLTCON MACH(3)=1.4$')
    assert result['MACH'] == [.6, None, 1.4]


def test_reassigned_scalar_stays_scalar():
    assert _flight('$FLTCON NMACH=3., NMACH=2.$')['NMACH'] == 2.


def test_ex2_complete_flight_schedules_reach_state():
    parser = NamelistParser()
    case = parser.parse_file(Path(__file__).parent / 'fixtures' / 'ex2.inp')[0]
    state = parser.to_state_dict(case)
    assert state['flight_mach'] == [.6, .9, 1.4, 2.5]
    assert state['flight_alt'] == [0., 2000., 40000., 90000.]
    assert state['flight_alschd'] == [-6., -4., -2., 0., 2., 4., 8., 12., 16., 20., 24.]
    assert len(state['flight_alschd']) == state['flight_nalpha']
