;Author Rigoberto Leander Salgado Reyes <rlsalgado2006 @gmail.com>
;
;Copyright 2016 by Rigoberto Leander Salgado Reyes.
;
;This program is licensed to you under the terms of version 3 of the
;GNU Affero General Public License. This program is distributed WITHOUT
;ANY EXPRESS OR IMPLIED WARRANTY, INCLUDING THOSE OF NON-INFRINGEMENT,
;MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. Please refer to the
;AGPL (http:www.gnu.org/licenses/agpl-3.0.txt) for more details.

(ns cuckoo-search-clj.cuckoo-search
  (:require [cuckoo-search-clj.operators :refer :all]))

(defn create-eggs
  "Function to create initials solutions
   nd   ->   number of dimensions.
   eggs ->   number of nests (or different solutions).
   gen  ->   function to generate random dimension.
   fun  ->   fitness function.
  "
  [nd eggs gen fun]
  (->> (fn [] (let [egg {:solution (->> gen (repeatedly nd) vec)}] (assoc egg :fitness (fun egg))))
    (repeatedly eggs)
    vec))

(defn create-cuckoo-search
  "Function to create initial state of the search.
   nd   ->   number of dimensions.
   eggs ->   number of nests (or different solutions).
   gen  ->   function to generate random dimension.
   fun  ->   fitness function.
   stop ->   stop function.
   lb   ->   lower bounds.
   ub   ->   upper bounds.
   pa   ->   discovery rate of alien eggs/solutions.
  "
  [nd eggs gen fun stop lb ub pa]
  (->> {:nest (create-eggs nd eggs gen fun)
        :new-nest []
        :fun fun
        :lb lb
        :ub ub
        :pa pa
        :perm1 (->> eggs (range 0) vec shuffle)
        :perm2 (->> eggs (range 0) vec shuffle)
        :stop stop
        :niter 0}
    (#(assoc % :best-nest (get-in % [:nest 0])))
    check-best-nest))

(defn search
  "search method
  [cs]        ->  cuckoo search initial state.
  [cs & ops]  ->  cs state and operators to process it.
  "
  ([cs] (search cs get-cuckoos best-nest empty-nest best-nest))
  ([{:keys [stop best-nest niter], :as cs*} & ops]
    ;    (prn (:fitness best-nest))
    (if (stop best-nest)
      {:best-nest best-nest :niter niter}
      (recur (update (reduce #(%2 %) cs* ops) :niter inc) ops))))

(defn random
  "Random value between lb and ub
  lb -> lower bounds
  up -> upper bounds
  "
  [lb ub]
  (+ (* (rand) (- ub lb)) lb))
