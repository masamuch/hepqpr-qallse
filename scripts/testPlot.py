from hepqpr.qallse import *
from hepqpr.qallse.plotting import * 


import time
import pickle

# import the method
from hepqpr.qallse.dsmaker import create_dataset

modelName = "D0"
#modelName = "Mp"
#modelName = "Doublet"

# 5e-3 : 167 MeV
# 8e-4 : 1.04 GeV
varDensity = {
    (modelName, 0.01, 5e-3),
    (modelName, 0.02, 5e-3),
    (modelName, 0.03, 5e-3),
    (modelName, 0.05, 5e-3),
    (modelName, 0.07, 5e-3),
    (modelName, 0.10, 5e-3),
    (modelName, 0.12, 5e-3),
    (modelName, 0.15, 5e-3),
    (modelName, 0.17, 5e-3),
    (modelName, 0.20, 5e-3),
    (modelName, 0.25, 5e-3),

    (modelName, 0.20, 8e-4),
    (modelName, 0.20, 1e-3),
    (modelName, 0.20, 2e-3),
    (modelName, 0.20, 3e-3),
}
var = {}
for v in varDensity:
    k = v[0]+f"_{v[1]*100:03.0f}perc_{v[2]:1.1e}"
    var[k] = v

print (var)

picklename = ".tmp.testPlot."+modelName+".pickle"
try:
    with open(picklename,'rb') as f:
        results = pickle.load(f)
except:
    print ("No picklet files.")
    results = {}

for k,v in var.items():
    Density = v[1]
    ModelName = v[0]
    ptThr = v[2]
    if k in results:
        results[k]["ptThr"] = ptThr
        results[k]["ModelName"] = ModelName
        continue

    results[k] = {}
    results[k]["density"] = Density
    results[k]["ptThr"] = ptThr
    results[k]["ModelName"] = ModelName
    # dataset creation options
    ds_options = dict(
        # output directory: output_path+prefix
        output_path='/tmp',
        prefix='ds_'+k,
        # size
        density = Density,
        phi_bounds = (0.15, 1.05),
        # important: no pt cut 
        high_pt_cut = 0,
    )

    # generate the dataset
    start_time = time.process_time()
    meta, path = create_dataset(**ds_options)
    exec_time = time.process_time() - start_time
    results[k]['TReadingHits'] = exec_time
    results[k]['meta']=meta

    from hepqpr.qallse.seeding import generate_doublets, SeedingConfig
    # generate the doublets: the important part is the config_cls !
    start_time = time.process_time()
    doublets = generate_doublets(hits_path=path+'-hits.csv', config_cls=SeedingConfig)
    exec_time = time.process_time() - start_time
    results[k]['TInitialDoubletBuilding'] = exec_time
    print('number of doublets = ', len(doublets))
    results[k]['Ndoublets'] = len(doublets)

    # optional: save the doublets to file
    doublets.to_csv(path+'-doublets.csv')


    if modelName == "D0":
        # Example: use QallseD0 with low pt
        # First, import the D0 config and the low pt base config
        from hepqpr.qallse.qallse_d0 import D0Config
        from hepqpr.qallse.qallse import Config

        # Just merge them: the second argument will override the first
        new_config = merge_dicts(D0Config().as_dict(), Config().as_dict())
        new_config["tplet_max_curv"] = ptThr

        # Now, just create the model using the new configuration
        dw = DataWrapper.from_path(path + '-hits.csv')
        model = QallseD0(dw, **new_config)
    elif modelName == "Mp":
        # Example: use QallseMp with low pt
        # First, import the Mp config and the low pt base config
        from hepqpr.qallse.qallse_mp import MpConfig
        from hepqpr.qallse.qallse import Config

        # Just merge them: the second argument will override the first
        new_config = merge_dicts(MpConfig().as_dict(), Config().as_dict())
        new_config["tplet_max_curv"] = ptThr

        # Now, just create the model using the new configuration
        dw = DataWrapper.from_path(path + '-hits.csv')
        model = QallseMp(dw, **new_config)
    elif modelName == "Nominal":
        # Example: use Qallse with low pt
        # First, import the 1GeV config and the low pt base config
        from hepqpr.qallse.qallse import Config, Config1GeV

        # Just merge them: the second argument will override the first
        new_config = merge_dicts(Config1GeV().as_dict(), Config().as_dict())
        new_config["tplet_max_curv"] = ptThr

        # Now, just create the model using the new configuration
        dw = DataWrapper.from_path(path + '-hits.csv')
        model = Qallse1GeV(dw, **new_config)
    elif modelName == "Doublet":
        # Example: use QallseDoublet with low pt
        # First, import the Doublet config and the low pt base config
        from hepqpr.qallse.qallse_doublet import DoubletConfig
        from hepqpr.qallse.qallse import Config

        # Just merge them: the second argument will override the first
        new_config = merge_dicts(DoubletConfig().as_dict(), Config().as_dict())
        new_config["tplet_max_curv"] = ptThr

        # Now, just create the model using the new configuration
        dw = DataWrapper.from_path(path + '-hits.csv')
        model = QallseDoublet(dw, **new_config)

    # generate the qubo as usual
    start_time = time.process_time()
    model.build_model(doublets)
    exec_time = time.process_time() - start_time
    print(f'Time of model building = {exec_time:.2f}s.')
    results[k]['TModelBuilding'] = exec_time
    Q = model.to_qubo()
    results[k]['QuboSize'] = len(Q)

    from hepqpr.qallse.cli.func import *
    start_time = time.process_time()
    response = solve_neal(Q)
    exec_time = time.process_time() - start_time
    print(f'Time of neal = {exec_time:.2f}s.')
    results[k]['TNeal'] = exec_time
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

import matplotlib.pyplot as plt
import numpy as np


variables = {
    "Density":[],
    "Nhits":[],
    "Ntracks":[],
    "Ndoublets":[],
    "TReadingHits":[],
    "TInitialDoubletBuilding":[],
    "TModelBuilding":[],
    "TNeal":[],
    "diffEnergy":[],
    "OccFrac":[],
    "Precision":[],
    "Recall":[],
    "Missing":[],
    "TrackmlScore":[],
    "QuboSize":[]

}

for k, v in sorted(results.items(), key=lambda x: x[1]["density"]):
    #print(v['meta'])
    #if "5.0e-03" not in k:
    if "020perc" not in k:
        continue
    variables["Density"]                 .append(v["density"]*100)
    variables["Nhits"]                   .append(v["meta"]["num_hits"])
    variables["Ntracks"]                 .append(v["meta"]["num_tracks"])
    variables["Ndoublets"]               .append(v["Ndoublets"])
    variables["TReadingHits"]            .append(v["TReadingHits"])
    variables["TInitialDoubletBuilding"] .append(v["TInitialDoubletBuilding"])
    variables["TModelBuilding"]          .append(v["TModelBuilding"])
    variables["TNeal"]                   .append(v["TNeal"])
    variables["diffEnergy"]              .append(v["obsEnergy"]-v['idealEnergy'])
    variables["OccFrac"]                 .append(v["bestOcc"]/float(v["OccSum"]))
    variables["Precision"]               .append(v["precision"]*100)
    variables["Recall"]                  .append(v["recall"]*100)
    variables["Missing"]                 .append(v["missing"])
    variables["TrackmlScore"]            .append(v["trackmlScore"])
    variables["QuboSize"]                .append(v["QuboSize"])

plotVar = {
    ("Ndoublets", "TModelBuilding"),
    ("Density"  , "TModelBuilding"),
    ("Density"  , "TReadingHits"),
    ("Density"  , "TInitialDoubletBuilding"),
    ("Density"  , "TNeal"),
    ("QuboSize" , "TNeal"),
    ("Density"  , "Precision"),
    ("Density"  , "Recall"),
    ("Density"  , "Missing"),
    ("Density"  , "diffEnergy"),
    ("Density"  , "Nhits"),
    ("Density"  , "Ntracks"),
    ("Density"  , "Ndoublets"),
    ("Density"  , "TrackmlScore")
}
legendStr = {
    "Density"       :r"Density [%]",
    "Nhits"         :r"N of hits",
    "Ntracks"       :r"N of tracks",
    "Ndoublets"     :r"N of doublets",
    "TReadingHits"  :r"Time of reading hits",
    "TInitialDoubletBuilding":r"Time of initial doublet building",
    "TModelBuilding":r"Time of model building",
    "TNeal"         :r"Time of SimulatedAnnealingSampler",
    "diffEnergy"    :r"Observed energy - ideal energy",
    "OccFrac"       :r"Best sample occurrence / sum of occurences",
    "Precision"     :r"Precision (Purity) [%]",
    "Recall"        :r"Recall (Efficienty) [%]",
    "Missing"       :r"N of missing doublets",
    "QuboSize"      :r"N of QUBO",
    "TrackmlScore"  :r"TrackML Score"
}
for v in plotVar:
    #print(v[0], v[1])
    #print(variables[v[1]])
    #plt.plot(variables[v[0]], variables[v[1]],marker='ro',label="data")
    plt.plot(variables[v[0]], variables[v[1]], 'ro', label="data")
    #plt.plot(variables[v[0]], np.poly1d(np.polyfit(variables[v[0]], variables[v[1]], 2))(variables[v[0]]), label='poly2')
    #plt.plot(variables[v[0]], np.poly1d(np.polyfit(variables[v[0]], variables[v[1]], 3))(variables[v[0]]), label='poly3')
    #plt.plot(variables[v[0]], np.poly1d(np.polyfit(variables[v[0]], variables[v[1]], 4))(variables[v[0]]), label='poly4')
    plt.xlabel(legendStr[v[0]])
    plt.ylabel(legendStr[v[1]])
    plt.legend()
    plt.savefig("fig_"+v[0]+"_"+v[1]+".png")
    plt.close()
    #plt.show()
