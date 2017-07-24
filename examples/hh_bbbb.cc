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
#include <momemta/Unused.h>

#include <TH1D.h>

#include <chrono>

using namespace std::chrono;
using namespace momemta;

int main(int argc, char** argv) {

    UNUSED(argc);
    UNUSED(argv);

    logging::set_level(logging::level::debug);

    ConfigurationReader configuration("../examples/hh_bbbb.lua");

    MoMEMta weight(configuration.freeze());

    // bjet
    Particle bjet1 { "bjet1", LorentzVector(15.68735, -189.8133, -210.1001, 288.65808), 5 };
    // anti-bjet
    Particle bjet2 { "bjet2", LorentzVector(-81.89723, -59.53946, -75.95126, 127.12102), -5 };
    // bjet
    Particle bjet3 { "bjet3", LorentzVector(84.728727, 127.45076, -243.5151, 289.04342), 5 };
    // anti-bjet
    Particle bjet4 { "bjet4", LorentzVector(-18.90508, 86.760458, -76.85704, 117.55769), -5 };

    auto start_time = system_clock::now();
    std::vector<std::pair<double, double>> weights = weight.computeWeights({bjet1, bjet2, bjet3, bjet4});
    auto end_time = system_clock::now();

    LOG(debug) << "Result:";
    for (const auto& r: weights) {
        LOG(debug) << r.first << " +- " << r.second;
    }

    LOG(info) << "Weight computed in " << std::chrono::duration_cast<milliseconds>(end_time - start_time).count() << "ms";


    return 0;
}
