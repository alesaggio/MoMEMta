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

#include <momemta/Logging.h>
#include <momemta/ParameterSet.h>
#include <momemta/Module.h>
#include <momemta/Types.h>
#include <momemta/Utils.h>

#include <TMath.h>

/** \brief \f$\require{cancel}\f$ Final (main) Block B, describing \f$q_1 q_2 \to X + s_{12} (\to \cancel{p_1} p_2)\f$
 *
 * \f$q_1\f$ and \f$q_2\f$ are Bjorken fractions, and \f$s_{12}\f$ is a particle decaying 
 * into \f$p_1\f$ (invisible particle) and \f$p_2\f$ (visible particle).
 *
 * This Block addresses the change of variables needed to pass from the standard phase-space
 * parametrization to the \f$\frac{1}{4\pi E_1} ds_{12} \times J\f$ parametrization.
 * 
 * The integration is performed over \f$s_{12}\f$ with \f$p_2\f$ as input. Per integration point, 
 * the LorentzVector of the invisible particle, \f$p_1\f$, is computed based on the following set 
 * of equations:   
 *
 * - \f$s_{12} = (p_1 + p_2)^2\f$
 * - Conservation of momentum (with \f$\vec{p}_T^{tot}\f$ the total transverse momentum of visible particles):
 *  - \f$p_{1x} = - p_{Tx}^{tot}\f$
 *  - \f$p_{1y} = - p_{Ty}^{tot}\f$
 * - \f$p_1^2 = m_1^2 = 0\f$ (FIXME)
 *
 * Up to two \f$p_1\f$ solutions are possible.
 *
 * ### Integration dimension
 *
 * This module adds **0** dimension to the integration.
 *
 * ### Global parameters
 *
 *   | Name | Type | %Description |
 *   |------|------|--------------|
 *   | `energy` | double | Collision energy. |
 *
 * ### Parameters
 *
 *   | Name | Type | %Description |
 *   |------|------|--------------|
 *   | `pT_is_met` | bool, default false | Fix \f$\vec{p}_{T}^{tot} = -\vec{\cancel{E_T}}\f$ or \f$\vec{p}_{T}^{tot} = \sum_{i \in \text{ vis}} \vec{p}_i\f$ |
 *
 * ### Inputs
 *
 *   | Name | Type | %Description |
 *   |------|------|--------------|
 *   | `s12` | double | Invariant mass of the particle decaying into the missing particle (\f$p_1\f$) and the visible particle, \f$p_2\f$. Typically coming from a BreitWignerGenerator module.
 *   | `inputs` | vector(LorentzVector) | LorentzVector of all the experimentally reconstructed particles. In this Block there is only one visible particle used explicitly, \f$p_2\f$, but there can be other visible objects in the the event, taken into account when computing \f$\vec{p}_{T}^{tot}\f$.
 *   | `met` | LorentzVector, default `input::met` | LorentzVector of the MET |
 *
 * ### Outputs
 *
 *   | Name | Type | %Description |
 *   |------|------|--------------|
 *   | `invisibles` | vector(vector(LorentzVector)) | LorentzVector of the invisible particles. In this Block \f$p_1\f$. One value per solution.
 *   | `jacobians` | vector(double) | Jacobian of the performed change of variables, leading to an integration on \f$ds_{12}\f$. One jacobian per solution.
 *
 * \warning This block is **not** validated for the moment. The output is maybe correct, maybe not. Use with caution.
 *
 * \ingroup modules
 */

class BlockA: public Module {
    public:
  
        BlockA(PoolPtr pool, const ParameterSet& parameters): Module(pool, parameters.getModuleName()) {

            sqrt_s = parameters.globalParameters().get<double>("energy");

            m_ps_point1 = parameters.get<InputTag>("theta1");
            m_ps_point1.resolve(pool);
            m_ps_point2 = parameters.get<InputTag>("phi1");
            m_ps_point2.resolve(pool);
            m_ps_point3 = parameters.get<InputTag>("theta2");
            m_ps_point3.resolve(pool);
            m_ps_point4 = parameters.get<InputTag>("phi2");
            m_ps_point4.resolve(pool);
                        
            m_particle_tags = parameters.get<std::vector<InputTag>>("inputs");
            for (auto& t: m_particle_tags)
              t.resolve(pool);
        }; 
  
        virtual void work() override {

            invisibles->clear();
            jacobians->clear();

            invisibles->push_back({});
            jacobians->push_back(computeJacobian());
       }

        // The extra dimensions not present in the input
        virtual size_t dimensions() const override {
            return 4;
        }
 
        double computeJacobian() {

            const LorentzVector& p3 = m_particle_tags[0].get<LorentzVector>();
            const LorentzVector& p4 = m_particle_tags[1].get<LorentzVector>();

            LorentzVector pb = p3+p4;
            for (size_t i = 2; i < m_particle_tags.size(); i++) {
                pb += m_particle_tags[i].get<LorentzVector>();
            }

            const double pbx = pb.Px();
            const double pby = pb.Py();

            double theta1 = m_ps_point1.get<double>()*M_PI;
            double phi1 = m_ps_point2.get<double>()*2*M_PI;
            double theta2 = m_ps_point3.get<double>()*M_PI;
            double phi2 = m_ps_point4.get<double>()*2*M_PI;

            //pT = p1+p2+pb. Balance it with the following system:
            //p1x+p2x = -pbx
            //p1y+p2y = -pby,
            //where: p1x=modp1*sin(theta1)*cos(phi1), p1y=modp1*sin(theta1)*sin(phi1),
            //       p2x=modp2*sin(theta2)*cos(phi2), p2y=modp2*sin(theta2)*sin(phi2)
            //Get modp1, modp2 as solutions of this system

            const double modp1 = -pbx/(std::sin(theta1)*std::cos(phi1)) - ((std::sin(theta2)*std::cos(phi2))*(pbx*std::tan(phi1)-pby))/(std::sin(theta1)*std::sin(theta2)*std::sin(phi2-phi1));
            const double modp2 = (pbx*std::sin(phi1)-pby*std::cos(phi1))/(std::sin(theta2)*std::sin(phi2-phi1));
   
            double inv_jac = 8*SQ(M_PI*sqrt_s)*(std::cos(phi1)*std::sin(phi2)-std::sin(phi1)*std::cos(phi2));
            
            return (SQ(modp1)*SQ(modp2)/(modp1*modp2)) / std::abs(inv_jac);    //WARNING: I'm putting E1=|p1|, E2=|p2|, we sholud consider E1=sqrt(|p|^2+m^2)...
        }

    private:
        double sqrt_s;

        std::vector<InputTag> m_particle_tags;

        InputTag m_ps_point1;
        InputTag m_ps_point2;
        InputTag m_ps_point3;
        InputTag m_ps_point4;

        std::shared_ptr<std::vector<std::vector<LorentzVector>>> invisibles = produce<std::vector<std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>>>>>("invisibles");
        std::shared_ptr<std::vector<double>> jacobians = produce<std::vector<double>>("jacobians");
};
REGISTER_MODULE(BlockA);
