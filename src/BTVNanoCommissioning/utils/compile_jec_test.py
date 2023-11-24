from BTVNanoCommissioning.utils.jet_factory import jet_factories

if __name__ == "__main__":
    import sys
    import gzip

    # jme stuff not pickleable in coffea
    import cloudpickle

    campaign = sys.argv[-2]
    jet_factory = jet_factories(campaign)
    print(campaign)

    with gzip.open(
        f"src/BTVNanoCommissioning/data/JME/{campaign}/{sys.argv[-1]}.pkl.gz", "wb"
    ) as fout:
        cloudpickle.dump(
            {
                "jet_factory": jet_factories(campaign),
                "met_factory": CorrectedMETFactory(jec_name_map),
            },
            fout,
        )
