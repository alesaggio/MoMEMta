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
#include <momemta/Logging.h>
#include <momemta/MoMEMta.h>
#include <momemta/Utils.h>

#include <TH1D.h>

#include <chrono>

using namespace std::chrono;

int main(int argc, char** argv) {

    UNUSED(argc);
    UNUSED(argv);

    logging::set_level(boost::log::trivial::debug);

    ConfigurationReader configuration("../examples/HH_bbbb.lua");

    MoMEMta weight(configuration.freeze());

    LorentzVector p1(16.171895980835, -13.7919054031372, -3.42997527122497, 21.5293197631836);
    LorentzVector p2(-18.9018573760986, 10.0896110534668, -0.602926552295686, 21.4346446990967);
    // b-quark
    LorentzVector p3(-55.7908325195313, -111.59294128418, -122.144721984863, 174.66259765625);
    // Anti b-quark
    LorentzVector p4(71.3899612426758, 96.0094833374023, -77.2513122558594, 142.492813110352);

    auto start_time = system_clock::now();
    std::vector<std::pair<double, double>> weights = weight.computeWeights({p1, p2, p3, p4});
    auto end_time = system_clock::now();

    LOG(debug) << "Result:";
    for (const auto& r: weights) {
        LOG(debug) << r.first << " +- " << r.second;
    }

    LOG(info) << "Weight computed in " << std::chrono::duration_cast<milliseconds>(end_time - start_time).count() << "ms";


    return 0;
}
