(defproject clj-mzml "0.1.1"
  :description "FIXME: write description"
  :url "http://example.com/FIXME"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [[org.clojure/clojure "1.8.0"]
                 [org.clojure/data.xml "0.0.8"]
                 [org.clojure/data.zip "0.1.2"]
                 [me.raynes/fs "1.4.6"]
                 [commons-codec/commons-codec "1.10"]
                 [biodb "0.2.3"]]
  :repl-options {:init (set! *print-length* 100)})
