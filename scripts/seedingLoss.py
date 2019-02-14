from hepqpr.qallse import *
from hepqpr.qallse.plotting import * 
from hepqpr.qallse.cli.func import time_this


import time
import pickle

# import the method
from hepqpr.qallse.dsmaker import create_dataset

from hepqpr.qallse.seeding import SeedingConfig
class nPhi53SeedingConfig(SeedingConfig):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.nPhiSlices = 53
        self._compute_derived_attrs()

class nPhi75SeedingConfig(SeedingConfig):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.nPhiSlices = 75
        self._compute_derived_attrs()

class nPhi100SeedingConfig(SeedingConfig):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.nPhiSlices = 100
        self._compute_derived_attrs()

class nPhi25SeedingConfig(SeedingConfig):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.nPhiSlices = 25
        self._compute_derived_attrs()

class nPhi10SeedingConfig(SeedingConfig):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.nPhiSlices = 10
        self._compute_derived_attrs()

class nPhi6SeedingConfig(SeedingConfig):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.nPhiSlices = 6
        self._compute_derived_attrs()

class nPhi1SeedingConfig(SeedingConfig):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.nPhiSlices = 1
        self._compute_derived_attrs()

modelName = "D0"
#modelName = "Mp"
#modelName = "Doublet"

maxTry=4

# 5e-3 : 167 MeV
# 8e-4 : 1.04 GeV
varDensity = [
    (modelName, 0.05, 5e-3, 53, maxTry),
    (modelName, 0.05, 5e-3, 100, maxTry),
    (modelName, 0.05, 5e-3, 75, maxTry),
    (modelName, 0.05, 5e-3, 25, maxTry),
    (modelName, 0.05, 5e-3, 10, maxTry),
    (modelName, 0.05, 5e-3, 6, maxTry),
    #(modelName, 0.05, 5e-3, 1, maxTry),
]

picklename = ".tmp.seedingLoss.pickle"
try:
    with open(picklename,'rb') as f:
        results = pickle.load(f)
except:
    print ("No pickle files.")
    results = {}

for v in varDensity:
    nTry = v[4]
    for iTry in range(nTry):
        k = (v[0], v[1], v[2], v[3], iTry)
        print (k)
        ModelName = k[0]
        Density   = k[1]
        ptThr     = k[2]
        nPhi      = k[3]
        if k in results:
            continue

        results[k] = {}
        results[k]["density"] = Density
        results[k]["ptThr"] = ptThr
        results[k]["ModelName"] = ModelName
        # dataset creation options
        ds_options = dict(
            # output directory: output_path+prefix
            output_path='/tmp',
            #prefix='ds_'+k,
            #prefix=prefix,
            # size
            density = Density,
            #phi_bounds = (0.15, 1.05),
            # important: no pt cut 
            high_pt_cut = 0,
        )

        prefix = f'ez-{Density}'
        if ds_options["high_pt_cut"] > 0:
            prefix += f'_hpt-{ds_options["high_pt_cut"]}'
        else:
            prefix += '_baby'
        prefix += f'_{iTry}'
        prefix += f'_nPhi{nPhi}'
        prefix += f'_noPhiCut'
        ds_options["prefix"] = prefix

        # generate the dataset
        import os
        path = os.path.join(ds_options['output_path'], prefix, "event000001000")
        if os.path.exists(path + "-hits.csv"):
            import json
            with open(path + "-meta.json") as f:
                meta = json.load(f)
            with open(path+"-metaHits.pickle", 'rb') as f:
                time_info= pickle.load(f)
        else:
            with time_this() as time_info:
                meta, path = create_dataset(**ds_options)
            with open(os.path.join(path+"-metaHits.pickle"), 'wb') as f:
                pickle.dump(time_info, f)

        results[k]['TReadingHits'] = time_info[1]
        results[k]['meta']=meta

        from hepqpr.qallse.seeding import generate_doublets, SeedingConfig
        if nPhi == 53  : seedingConfig = nPhi53SeedingConfig
        if nPhi == 100 : seedingConfig = nPhi100SeedingConfig
        if nPhi == 75  : seedingConfig = nPhi75SeedingConfig
        if nPhi == 25  : seedingConfig = nPhi25SeedingConfig
        if nPhi == 10  : seedingConfig = nPhi10SeedingConfig
        if nPhi == 6   : seedingConfig = nPhi6SeedingConfig
        if nPhi == 1   : seedingConfig = nPhi1SeedingConfig
        # generate the doublets: the important part is the config_cls !
        if os.path.exists(path + "-doublets.csv"):
            doublets = pd.read_csv(path + "-doublets.csv", index_col=0)
            results[k]['TInitialDoubletBuilding'] = time_info[1]
            with open(path+"-metaDoublets.pickle", 'rb') as f:
                time_info= pickle.load(f)
        else:
            with time_this() as time_info:
                doublets = generate_doublets(hits_path=path+'-hits.csv', config_cls=seedingConfig)
                doublets.to_csv(path+'-doublets.csv')
            with open(os.path.join(path+"-metaDoublets.pickle"), 'wb') as f:
                pickle.dump(time_info, f)
        results[k]['TInitialDoubletBuilding'] = time_info[1]
        print('number of doublets = ', len(doublets))
        results[k]['Ndoublets'] = len(doublets)



        from hepqpr.qallse.qallse import Config
        config = Config()
        config.tplet_max_curv = ptThr
        dw = DataWrapper.from_path(path + '-hits.csv')
        if modelName == "D0":
            from hepqpr.qallse.qallse_d0 import D0Config
            new_config = merge_dicts(D0Config().as_dict(), config.as_dict())
            model = QallseD0(dw, **new_config)
        elif modelName == "Mp":
            from hepqpr.qallse.qallse_mp import MpConfig
            new_config = merge_dicts(MpConfig().as_dict(), config.as_dict())
            model = QallseMp(dw, **new_config)
        elif modelName == "Nominal":
            from hepqpr.qallse.qallse import Config1GeV
            new_config = merge_dicts(Config1GeV().as_dict(), config.as_dict())
            model = Qallse1GeV(dw, **new_config)
        elif modelName == "Doublet":
            from hepqpr.qallse.qallse_doublet import DoubletConfig
            new_config = merge_dicts(DoubletConfig().as_dict(), config.as_dict())
            model = QallseDoublet(dw, **new_config)

        p, r, ms = model.dataw.compute_score(doublets)
        results[k]['precision_initDoublet'] = p
        results[k]['recall_initDoublet'] = r
        results[k]['missing_initDoublet'] = len(ms)

        # generate the qubo as usual
        with time_this() as time_info:
            model.build_model(doublets)
        print(f'Time of model building = {time_info[1]:.2f}s.')
        results[k]['TModelBuilding'] = time_info[1]

        with time_this() as time_info:
            Q = model.to_qubo()
        print(f'Time of qubo building = {time_info[1]:.2f}s.')
        results[k]['TQuboBuilding'] = time_info[1]
        results[k]['QuboSize'] = len(Q)

        from hepqpr.qallse.cli.func import *
        with time_this() as time_info:
            response = solve_neal(Q)
        print(f'Time of neal = {time_info[1]:.2f}s.')
        results[k]['TNeal'] = time_info[1]
        final_doublets, final_tracks = process_response(response)

        en0 = 0 if Q is None else dw.compute_energy(Q)
        en = response.record.energy[0]
        results[k]['obsEnergy'] = en
        results[k]['idealEnergy'] = en0
        occs = response.record.num_occurrences
        results[k]['bestOcc'] = occs[0]
        results[k]['OccSum'] = occs.sum()

        p, r, ms = dw.compute_score(final_doublets)
        results[k]['precision'] = p
        results[k]['recall'] = r
        results[k]['missing'] = len(ms)
        trackml_score = dw.compute_trackml_score(final_tracks)
        results[k]['trackmlScore'] = trackml_score


with open(picklename, 'wb') as f:
    pickle.dump(results, f)

#print(results)
