'''
created: oct 2020
last modified: oct 2020
owner: bill

bio: test program for illustris_api.get()
note: copied from https://www.illustris-project.org/data/docs/api/
'''

import requests

baseUrl = 'http://www.illustris-project.org/api/'
headers = {"api-key":"e5d5b9ee9a836b33d764bb21d88d178c"}

def get(path=baseUrl, params=None, prefix=''):
    # make HTTP GET request to path
    r = requests.get(path, params=params, headers=headers)

    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        return r.json() # parse json responses automatically
    
    if 'content-disposition' in r.headers:
        filename = prefix + r.headers['content-disposition'].split("filename=")[1]
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string

    return r