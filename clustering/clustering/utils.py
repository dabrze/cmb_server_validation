# coding: utf-8
__author__ = 'Dariusz Brzezinski'

import psycopg2
import logging
import sys
import os
import requests
import urllib
import gzip
import shutil

from sqlalchemy import create_engine

from abc import ABC, abstractmethod


class DbScript(ABC):
    def __init__(self, name, database_name='pdbmonomers', port=55543, tables_to_drop=[], logging_level=logging.INFO,
                 log_file=None):
        self.script = ""
        self.name = name
        self.database_name = database_name
        self.port = port
        self.conn = None
        self.cursor = None
        self.tables_to_drop = tables_to_drop
        self._setup_logging(log_file, logging_level)

    def _setup_logging(self, log_file, logging_level):
        self.logger = logging.getLogger("DbScript")
        self.logger.setLevel(logging_level)

        if log_file is not None:
            handler = logging.FileHandler('logfile.log')
        else:
            handler = logging.StreamHandler(sys.stdout)

        formatter = logging.Formatter('%(asctime)s:%(levelname)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        handler.setFormatter(formatter)
        self.logger.addHandler(handler)

    def setup_connection(self):
        self.conn = psycopg2.connect(database=self.database_name, port=self.port)
        self.cursor = self.conn.cursor()

    def drop_tables(self):
        drop_script = "\n".join(["drop table if exists {0};".format(t) for t in self.tables_to_drop])
        if drop_script != "":
            self.logger.info("Droping tables: {0}".format(", ".join(self.tables_to_drop)))
            self.cursor.execute(drop_script)

    def set_functions(self):
        self.cursor.execute("""
            CREATE OR REPLACE FUNCTION array_sort (ANYARRAY)
             RETURNS ANYARRAY LANGUAGE SQL
             AS $$SELECT ARRAY(SELECT unnest($1) ORDER BY 1)$$;
            
            create extension if not exists rdkit;
        """)

    def run_commands(self, commands):
        self.logger.debug(commands)
        self.cursor.execute(commands)
        self.conn.commit()

    def query_pdbj(self, sql, format, dest_file):
        pdbj_rest_url = "https://pdbj.org/rest/mine2_sql"
        params = {
            "q": sql,
            "format": format,
        }

        response = requests.get(pdbj_rest_url, params)
        response.raise_for_status()

        directory = os.path.dirname(dest_file)
        if not os.path.exists(directory):
            os.makedirs(directory)

        with open(dest_file, 'wb') as handle:
            for block in response.iter_content(2048):
                handle.write(block)

    def download_file(self, url, dest_file, unpack):
        directory = os.path.dirname(dest_file)
        if not os.path.exists(directory):
            os.makedirs(directory)

        if unpack:
            tmp_file = dest_file + ".gz"
            urllib.request.urlretrieve(url, tmp_file)

            with gzip.open(tmp_file, "rb") as gz_file:
                with open(dest_file, 'wb') as unpacked_file:
                    shutil.copyfileobj(gz_file, unpacked_file)
            os.remove(tmp_file)
        else:
            urllib.request.urlretrieve(url, dest_file + ".gz")

    def copy_file_to_db(self, file, table_name, columns=""):
        self.run_commands("""
            copy {table_name} {columns} from '{file}' with delimiter E'\t' CSV HEADER;
        """.format(table_name=table_name, columns=columns, file=file))

    def _get_conn(self):
        return self.conn

    def save_result_df(self, df, output_name, output_folder="data"):
        output_file = os.path.join(output_folder, output_name + ".tsv")
        engine = create_engine("postgresql+psycopg2://", creator=self._get_conn)

        self.logger.debug("Saving results to {0} and table {1}".format(output_file, output_name))
        df.to_csv(output_file, sep="\t")
        df.to_sql(output_name, engine, if_exists="replace")

    @abstractmethod
    def run(self):
        self.setup_connection()
        self.set_functions()
        self.drop_tables()

