#!/usr/bin/env python3

"""
Tests for count_multibase_positions in confindr_src.methods
"""

# Local imports
from confindr_src import methods


def test_count_multibase_positions_excludes_meta_keys():
    """
    Test that count_multibase_positions correctly counts the number of
    positions with multiple bases, excluding metadata keys.
    """
    # Simulate read_contig outputs
    multibase_dict_list = [
        {
            'geneA': {
                '_gene_stats': {
                    'combined_p': None,
                    'num_sig_positions': 0
                }
            }
        },
        {
            'geneB': {
                10: {'paired': {}},
                '_gene_stats': {'num_sig_positions': 0}
            }
        },
        {
            'geneC': {
                20: {'paired': {}},
                30: {'paired': {}},
                '_gene_stats': {'num_sig_positions': 2}
            }
        }
    ]

    total = methods.count_multibase_positions(
        multibase_dict_list=multibase_dict_list
    )
    # Expect 0 + 1 + 2 = 3
    assert total == 3
