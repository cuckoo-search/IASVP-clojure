;Author Rigoberto Leander Salgado Reyes <rlsalgado2006 @gmail.com>
;
;Copyright 2016 by Rigoberto Leander Salgado Reyes.
;
;This program is licensed to you under the terms of version 3 of the
;GNU Affero General Public License. This program is distributed WITHOUT
;ANY EXPRESS OR IMPLIED WARRANTY, INCLUDING THOSE OF NON-INFRINGEMENT,
;MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. Please refer to the
;AGPL (http:www.gnu.org/licenses/agpl-3.0.txt) for more details.

(ns cuckoo-search-clj.core
  (:gen-class)
  (:require [cuckoo-search-clj.cuckoo-search :refer :all]))

(defn d-dimensional-sphere
  []
  (search (create-cuckoo-search
            15 ; nd
            25 ; eggs
            (partial random -5 5) ; gen
            (fn [{solution :solution}] ; fun: d-dimensional sphere function sum_j=1^d (u_j-1)^2, with a minimum at (1,1, ...., 1)
              (reduce #(+ % (* (- %2 1) (- %2 1))) 0 solution))
            (fn [{fitness :fitness}] (<= fitness 1.0e-5)) ; stop function: when fitness(best-nest) < 1.0e-5
            -5 ; lb
            5 ; ub
            0.25 ; pa
            )))

(defn -main
  ""
  [& args]
  (println "Hello, World!"))



