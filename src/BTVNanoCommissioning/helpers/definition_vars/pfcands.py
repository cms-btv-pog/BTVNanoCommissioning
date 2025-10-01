def get_variables():
    vardict = {
        "PFCands_numberOFHits": {
            "displayname": "Number of hits",
            "ylabel_text": "Jets",
            "format_unit": "2f",
            "format_unit_digits": 2,
            "bins": 50,
            "manual_ranges": [0, 50],
            "inputVar_units": None,
        },
        "PFCands_numberOfPixelHits": {
            "displayname": "Number of pixel hits",
            "ylabel_text": "Jets",
            "format_unit": "2f",
            "format_unit_digits": 2,
            "bins": 10,
            "manual_ranges": [0, 10],
            "inputVar_units": None,
        },
        "PFCands_lostInnerHits": {
            "displayname": "Lost inner pixel hits",
            "ylabel_text": "Jets",
            "format_unit": "2f",
            "format_unit_digits": 2,
            "bins": 11,
            "manual_ranges": [-1, 10],
            "inputVar_units": None,
        },
    }

    return vardict
