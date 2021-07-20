# -*- coding: utf-8 -*-
"""
Author: Ying Huang

Date: 2020-03-17 23:14:25
Last Modified by: Ying Huang
Last Modified time: 2020-03-17 23:14:25

Description: 
Run Panther by using selenium
"""
import os
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.select import Select


PWD = os.path.split(os.path.realpath(__file__))[0]
Path_webdrivers = os.path.join(PWD, 'webdrivers')
Path_extension = os.path.join(PWD, 'extensions')

def wait(driver, time, **elements):
    """
    Wait for elements to be completely loaded.

    Required:
    driver: <object> object of Class 'webdriver.Chrome'.
    time: <int> wait time (seconds).

    Operation:
    elements: <dict> a dict for elements which keys are element names,
        values are element types in 'By' object.
    """
    print('<Link "{}", Wait for {}s each time>'.format(driver.current_url, time))
    ele_type_name = {
        By.ID: 'id',
        By.XPATH: 'xpath',
    }

    for i in range(10):
        try:
            wait = WebDriverWait(driver, time)
            break
        except KeyboardInterrupt:
            driver.quit()
        except:
            if i == 9:
                print('link reach {}, driver quit'.format(str(i)))
                driver.quit() 
                raise Exception('<ERR>: Can not link to "{}"'.format(driver.current_url))

            driver.refresh()
            print(
                '\t<Warning>: Failed to load website in "{}s", retry (refresh and wait) for "{}" times'.format(
                    str(time), str(i + 1))
            )

    for i in range(10):
        num_err_signal = 0

        print('\t<Checking elements {} time> ...'.format(str(i+1)))
        try:
            for j, (element, ele_type) in enumerate(elements.items()):
                try:
                    wait.until(EC.presence_of_element_located((ele_type, element)))
                except:
                    num_err_signal +=1
                    print(
                        '\t<Warning> : can not load element "{}" which type is "{}" in {}s'.format(
                            element, ele_type_name[ele_type], time
                        )
                    )
                    if j == len(elements.items()) - 1:
                        raise Exception(
                        '\t<Warning> : can not load element "{}" which type is "{}" in {}s'.format(
                            element, ele_type_name[ele_type], time
                        )
                    )

        except:
            if i == 9: # at last refresh time, can not locate element
                driver.quit()
                raise Exception('element check reach {}, driver quit'.format(str(i)))    
            else:
                driver.refresh()
                print('\t>Warning: failed to locate elements, retry for {} time'.format(str(i + 1)))
                continue
        if num_err_signal == 0:
            break

    return driver


def start_chrome(download_dir, if_headless=False, webdrivers='chromedriver.v80.exe', extensions='Ghelper_1.4.6.crx'):
    download_dir = os.path.abspath(download_dir)
    # set path to extensions and webdrivers
    Ghelper = os.path.join(Path_extension, extensions)
    print('> Using extension Ghelper: ', Ghelper)
    chromedriver = os.path.join(Path_webdrivers, webdrivers)
    print('> Using chromedriver: ', chromedriver)

    # set download place
    prefs = {
        'profile.default_content_settings.popups': 0, # shutdown pop-up windows
        'download.default_directory': download_dir, # set download directory
    }

    # add chrome options to webdriver
    chrome_options = webdriver.ChromeOptions()
    
    if if_headless is True:
        chrome_options.add_argument("--headless")
        chrome_options.add_argument('--disable-gpu')
        chrome_options.add_argument('--no-sandbox')

    chrome_options.add_extension(Ghelper)
    chrome_options.add_experimental_option('prefs', prefs) # add settings to browser

    # start chromedriver
    driver = webdriver.Chrome(
        executable_path=chromedriver,
        chrome_options=chrome_options,
    )

    return driver

def login_Ghelper(driver, url='chrome-extension://cieikaeocafmceoapfogpffaalkncpkc/login.html', email='hhhyyy1992117@qq.com', passwd='******'):

    driver.get(url)
    print('<Login Ghelper> ...')
    driver = wait(driver, 600, email=By.ID, password=By.ID, submit=By.ID) # wait 10 min
    
    # login Ghelper
    driver.find_element_by_id('email').send_keys(email)
    driver.find_element_by_id('password').send_keys(passwd)
    driver.find_element_by_id('submit').click()

    return driver

def start_panther_analysis(driver, geneids, odir, species='gal', annot='bp', url='http://pantherdb.org/'):
    species = species.lower()

    species_dict = {
        'gal': 'Gallus gallus',
        'hg': 'Homo sapiens',
        'cow': 'Bos taurus',
        'mus': 'Mus musculus',
    }

    annot = annot.lower()
    annot_dict = {
        'bp': 'fullgo_bp_comp',
        'cc': 'fullgo_cc_comp',
        'mf': 'fullgo_mf_comp',
        'pathway': 'pathway',
        'reactome': 'reactome',
    }

    print('<Start panther website with URL> :', url)
    driver.get(url)

    wait(
        driver, 10, # wait for 10s
        ** {
            'idField': By.ID, # find input region
            '//*[@id="mainBody"]/table[1]/tbody/tr[2]/td/div/form/table[1]/tbody/tr/td[2]/table/tbody/tr/td/table[3]/tbody/tr[2]/td[2]/table/tbody/tr[3]/td/input': By.XPATH, # select 'Statistical overrepresentation test'
            '//*[@id="mainBody"]/table[1]/tbody/tr[2]/td/div/form/table[2]/tbody/tr/td/a': By.XPATH, # 'submit' button
        }
    )

    driver.find_element_by_id('idField').send_keys(geneids)
    
    # select species
    wait(driver, 10, dataset=By.NAME)
    spec_select = Select(driver.find_element_by_name('dataset'))
    spec_select.deselect_all()
    spec_select.select_by_value(species_dict[species])

    print('> species: "{}", given "{}"'.format(spec_select.first_selected_option.text, species_dict[species]))
    
    # select analysis
    enrich_analysis = driver.find_element_by_xpath('//*[@id="mainBody"]/table[1]/tbody/tr[2]/td/div/form/table[1]/tbody/tr/td[2]/table/tbody/tr/td/table[3]/tbody/tr[2]/td[2]/table/tbody/tr[3]/td/input')
    enrich_analysis.click()
    wait(driver, 10, annotDataSet=By.ID)
    annot_select = Select(driver.find_element_by_id('annotDataSet')) 
    annot_select.select_by_value(annot_dict[annot])# select full go bp
    #print('<Select analysis is>: ', annot_select.first_selected_option.text)

    # submit
    driver.find_element_by_xpath('//*[@id="mainBody"]/table[1]/tbody/tr[2]/td/div/form/table[2]/tbody/tr/td/a').click()

    # select whole genome genes as reference genes
    wait(driver, 10, lists=By.NAME)
    ref_genes_select = Select(driver.find_element_by_name('lists'))
    ref_genes_select.select_by_value(species_dict[species])
    #print('<Select reference genes species is>: ', ref_genes_select.first_selected_option.text)

    #submit to run analysis
    wait(driver, 10, **{'//*[@id="mainBody"]/form/table/tbody/tr/td/table/tbody/tr[14]/td/a': By.XPATH})
    driver.find_element_by_xpath('//*[@id="mainBody"]/form/table/tbody/tr/td/table/tbody/tr[14]/td/a').click()

    # download result in JSON format
    wait(driver, 10, **{'/html/body/a[3]': By.XPATH})
    driver.find_element_by_xpath('/html/body/a[3]').click()

    return driver


def enrich_analysis(geneids, odir, species='gal', annot='bp', url='http://pantherdb.org/', if_headless=False, webdrivers='chromedriver.v80.exe', extensions='Ghelper_1.4.6.crx'):
    try:
        driver = start_chrome(odir, if_headless, webdrivers, extensions)
        driver = login_Ghelper(driver)
        driver = start_panther_analysis(driver, geneids, odir, species, annot, url)
    finally:
        driver.quit()
