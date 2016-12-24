(ns clj-mzml.core
  (:require [clojure.java.io :as io]
            [clojure.data.xml :refer [parse]]
            [clojure.zip :refer [node xml-zip]]
            [clojure.data.zip.xml :refer [xml-> attr xml1-> text attr=]]
            [biodb.core :as bdb])
  (:import [org.apache.commons.codec.binary Base64]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; utilities
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- parse-mzml
  [r]
  (let [x (parse r)]
    (if (= (:tag x) :indexedmzML)
      (first (:content x))
      x)))

(defn- get-parameter-tag
  [xml tag]
  (->> (take-while #(not (= (:tag %) :run))
                   (:content xml))
       (filter #(= (:tag %) tag))
       first
       xml-zip))

(defn- parse-cvs
  [z]
  (if (seq z)
    (->> (filter (fn [[k v]] (seq v))
                 {:cvparams (->> (map #(let [a (:attrs %)]
                                         (vector (:accession a) a))
                                      (xml-> z :cvParam node))
                                 (into {}))
                  :userparams (map :attrs (xml-> z :userParam node))
                  :references (->> (xml-> z :referenceableParamGroupRef (attr :ref)))})
         (into {}))))

(defn- parse-reference-group
  [zs]
  (->> (map #(vector (xml1-> % (attr :id))
                     (merge (-> (node %) :attrs (dissoc :id)) (parse-cvs %)))
            zs)
       (into {})))

(defn- components
  [z tag]
  (->> (xml-> z :componentList tag)
       (map #(merge {:order (xml1-> % (attr :order))} (parse-cvs %)))))

(defn- attr-cvs
  [z]
  (merge (-> z node :attrs) (parse-cvs z)))

(defn- inflate 
  [s]
  (let [h (new java.util.zip.InflaterInputStream
               (new java.io.ByteArrayInputStream s)
               (new java.util.zip.Inflater))]
    (try
      (loop [st h
             acc nil]
        (if (= (.available h) 0)
          (byte-array (map unchecked-byte (reverse acc)))
          (recur h (cons (.read h) acc))))
      (catch Exception e
        (println e))
      (finally
        (.close h)))))

(defn- get-buffer
  [cvs inflated]
  (let [btypes {"MS:1000519" #(.asIntBuffer inflated)
                "MS:1000521" #(.asFloatBuffer inflated)
                "MS:1000523" #(.asDoubleBuffer inflated)
                "MS:1000522" #(.asLongBuffer inflated)}
        accession (filter btypes cvs)]
    (cond (not (nil? (next accession)))
          (throw (Exception. "Multiple encoding cvParams."))
          (nil? (seq accession))
          (throw (Exception. "No supported encoding cvParam."))
          :else
          ((btypes (first accession))))))

(defn decode-binary
  [z]
  (let [decoded (.decode (Base64.) (xml1-> z :binary text))
        inflated (-> (if (xml1-> z :cvParam (attr= :accession "MS:1000576"))
                       decoded (inflate decoded))
                     java.nio.ByteBuffer/wrap
                     (.order java.nio.ByteOrder/LITTLE_ENDIAN))
        buffer (get-buffer (xml-> z :cvParam (attr :accession)) inflated)]
    {:array (loop [d buffer
                   acc nil]
              (if (not (.hasRemaining d))
                (vec (reverse acc))
                (recur d (cons (.get d) acc))))}))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; spectrum accessors
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn spectra-seq
  [reader]
  (let [xml (parse-mzml reader)]
    (->> (filter #(#{:run} (:tag %)) (:content xml)) first :content
         (filter #(#{:spectrumList} (:tag %))) first :content
         (filter #(#{:spectrum} (:tag %))))))

(defn index
  [s]
  "Returns the index of a spectrum."
  (-> s :attrs :index))

(defn id
  [s]
  "Returns the id of a spectrum."
  (-> s :attrs :id))

(defn- get-array
  [s acc]
  (->> (xml-> (xml-zip s) :binaryDataArrayList :binaryDataArray)
       (filter #(xml1-> % :cvParam (attr= :accession acc)))
       (map #(merge {:units (xml1-> % :cvParam (attr= :accession acc)
                                    (attr :unitName))}
                    (decode-binary %)))))

(defn mz-array
  [s]
  "Returns a list of maps containing the units and decoded mz arrays,
  keywords :units and :array respectively, from a spectrum."
  (get-array s "MS:1000514"))

(defn intensity-array
  [s]
  "Returns a list of maps containing the units and decoded intensity arrays,
  keywords :units and :array respectively, from a spectrum."
  (get-array s "MS:1000515"))

(defn intensity-array-units
  [s]
  "Returns a map of unit information for ")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; getting parameters, cv lists etc
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def parse-elements
  {:cvlist
   (fn [x] (->> (xml-> (get-parameter-tag x :cvList) :cv node) (map :attrs)))
   
   :file-content
   (fn [x] (->> (xml1-> (get-parameter-tag x :fileDescription) :fileContent)
               parse-cvs))

   :source-files
   (fn [x]
     (->> (xml-> (get-parameter-tag x :fileDescription)
                 :sourceFileList :sourceFile)
          parse-reference-group))

   :contacts
   (fn [x] (->> (xml-> (get-parameter-tag x :fileDescription) :contact)
               (map parse-cvs)))

   :referenceable
   (fn [x] (->> (xml-> (get-parameter-tag x :referenceableParamGroupList)
                      :referenceableParamGroup)
               parse-reference-group))

   :samples
   (fn [x] (->> (xml-> (get-parameter-tag x :sampleList) :sample)
               parse-reference-group))

   :software
   (fn [x] (->> (xml-> (get-parameter-tag x :softwareList) :software)
               parse-reference-group))

   :scan-settings
   (fn [x] (->> (xml-> (get-parameter-tag x :scanSettingsList) :scan-settings)
               (map #(vector (xml1-> % (attr :id))
                             (merge (-> (node %) :attrs (dissoc :id))
                                    (parse-cvs %)
                                    {:source-refs (xml-> % :sourceFileRefList
                                                         :sourceFileRef
                                                         (attr :ref))
                                     :targets (->> (xml-> % :targetList :target)
                                                   (map parse-cvs))})))))
   :instrument-config
   (fn [x] (->> (xml-> (get-parameter-tag x :instrumentConfigurationList)
                      :instrumentConfiguration)
               (map #(vector (xml1-> % (attr :id))
                             (merge (-> (node %) :attrs (dissoc :id))
                                    (parse-cvs %)
                                    {:sources (components % :source)
                                     :analyzers (components % :analyzer)
                                     :detectors (components % :detector)})))
               (into {})))

   :data-processing
   (fn [x] (->> (xml-> (get-parameter-tag x :dataProcessingList) :dataProcessing)
               (map #(vector (xml1-> % (attr :id))
                             (map (fn [x]
                                    (merge (-> x node :attrs)
                                           (parse-cvs x)))
                                  (xml-> % :processingMethod))))
               (into {})))

   :default-data-proc-ref
   (fn [x] (->> (filter #(#{:run} (:tag %)) (:content x))
               first :content first :attrs :defaultDataProcessingRef))

   :spectrum-count
   (fn [x] (Integer/parseInt (->> (filter #(#{:run} (:tag %)) (:content x))
                                 first :content first :attrs :count)))
   
   :run-parameter
   (fn [x]
     (let [re (->> (filter #(#{:run} (:tag %)) (:content x)) first)
           z (-> (assoc re :content
                        (take-while #(not (#{:spectrumList :chromatogramList}
                                           (:tag %)))
                                    (:content re)))
                 xml-zip)]
       (merge (-> z node :attrs) (parse-cvs z))))})

(defn mzml-parameters
  ([file] (mzml-parameters file false))
  ([file ks]
   (with-open [r (io/reader file)]
     (let [xml (parse-mzml r)]
       (->> (map (fn [[k v]] (vector k (v xml)))
                 (if ks 
                   (select-keys parse-elements ks)
                   parse-elements))
            (into {}))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; biodb compatability
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod bdb/table-spec :mzml
  [q]
  (vector [:accession :text "PRIMARY KEY"]
          [:number :integer "UNIQUE"]
          [:src :binary "NOT NULL"]))

(defmethod bdb/prep-sequences :mzml
  [q]
  (->> (:coll q)
       (map #(hash-map :accession (id %) :number (Integer/parseInt (index %))
                       :src (bdb/freeze %)))))

(defmethod bdb/restore-sequence :mzml
  [q]
  (bdb/thaw (:src (dissoc q :type))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; index mzmls
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn index-mzml-file
  [file]
  (let [db (bdb/db-spec {:dbname (str file ".sqlite")
                         :dbtype :sqlite})]
    (bdb/create-table! db :spectra :mzml)
    (with-open [r (io/reader file)]
      (bdb/insert-sequences! db :spectra :mzml (spectra-seq r)))
    db))

(defn get-spectra-by-index
  [db index]
  (bdb/query-sequences db ["select * from spectra where number=?" index] :mzml))
