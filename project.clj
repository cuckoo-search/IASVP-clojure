;Author Rigoberto Leander Salgado Reyes <rlsalgado2006 @gmail.com>
;
;Copyright 2016 by Rigoberto Leander Salgado Reyes.
;
;This program is licensed to you under the terms of version 3 of the
;GNU Affero General Public License. This program is distributed WITHOUT
;ANY EXPRESS OR IMPLIED WARRANTY, INCLUDING THOSE OF NON-INFRINGEMENT,
;MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. Please refer to the
;AGPL (http:www.gnu.org/licenses/agpl-3.0.txt) for more details.


(defproject cuckoo-search-clj "0.1.0"
  :description "Cuckoo Search Algorithm and
  Hybrid Cuckoo Search applied to the Inverse Additive Singular Values Problem(IASVP)"
  :url "https://www.researchgate.net/publication/268278164_Cuckoo_Search_hibrido_aplicado_al_Problema_Inverso_Aditivo_de_Valores_Singulares_Hybrid_Cuckoo_Search_applied_to_the_Inverse_Additive_Singular_Values_Problem"
  :license {:name "GNU Affero General Public License"
            :url "http:www.gnu.org/licenses/agpl-3.0.txt"}

  :dependencies [[org.clojure/clojure "1.8.0"]
                 [com.googlecode.matrix-toolkits-java/mtj "1.0.4"]]

  :main ^:skip-aot cuckoo-search-clj.core
  :target-path "target/%s"
  :profiles {:uberjar {:aot :all}})
