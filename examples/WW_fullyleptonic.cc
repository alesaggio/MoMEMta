/*
 *  MoMEMta: a modular implementation of the Matrix Element Method
 *  Copyright (C) 2016  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <momemta/ConfigurationReader.h>
#include <momemta/MoMEMta.h>
#include <momemta/Utils.h>

#include <TH1D.h>

#include <chrono>

using namespace std::chrono;

int main(int argc, char** argv) {

    UNUSED(argc);
    UNUSED(argv);

    spdlog::set_level(spdlog::level::trace);

    ConfigurationReader configuration("../examples/WW_fullyleptonic.lua");

    // Change top mass
    //configuration.getGlobalParameters().set("top_mass", 173.);

    MoMEMta weight(configuration.freeze());

    // Electron
    LorentzVector p3(16.171895980835, -13.7919054031372, -3.42997527122497, 21.5293197631836);
    // Muon
    LorentzVector p4(-18.9018573760986, 10.0896110534668, -0.602926552295686, 21.4346446990967);

    auto start_time = system_clock::now();
    std::vector<std::pair<double, double>> weights = weight.computeWeights({p3, p4});
    auto end_time = system_clock::now();

    LOG(debug) << "Result:";
    for (const auto& r: weights) {
        LOG(debug) << r.first << " +- " << r.second;
    }
    
    LOG(debug) << "Hist in pool: " << weight.getPool().exists({"dmem_WW", "hist"});
    std::shared_ptr<const TH1D> dmem = weight.getPool().get<TH1D>({"dmem_WW", "hist"});
    LOG(debug) << "DMEM integral: " << dmem->Integral();

    LOG(info) << "Weight computed in " << std::chrono::duration_cast<milliseconds>(end_time - start_time).count() << "ms";


    return 0;
}