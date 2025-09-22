""" This module contains routines to talk with NOMAD. They are only
slightly modified routines provided by the NOMAD guys.
"""
from ase2sprkkr.common.decorators import add_to_signature, cached_class_property
from functools import wraps


def with_token(func):

    @add_to_signature(func)
    @wraps(func)
    def func_with_token(self, token:str=None, *args, **kwargs):
        token = token or self.token
        if not token:
            raise ValueError("Nomad authentication token not set")
        func(*args, **kwargs, token=token)

    return func_with_token


class NomadApi():

    default_api_url='https://nomad-lab.eu/prod/v1/api/v1/'

    @cached_class_property
    def requests():
        """ Loading requests lasts for ages """
        import requests
        return requests

    def __init__(self, nomad_url=None):
        self.url = nomad_url or self.default_api_url
        self.token = None

    def get_authentication_token(self, username, password, expires=None):
        '''Get the token for accessing your NOMAD unpublished uploads remotely'''
        try:
            response = self.requests.get(
                self.url + 'auth/token', params=dict(username=username, password=password), timeout=10)
            token = response.json().get('access_token')
            if token:
                if expires:
                    response = self.requests.get(
                        self.url + 'auth/token', params=dict(token=token, username=username, password=password, expires=expires), timeout=10)
                    token = response.json().get('access_token')

                if token:
                    return token

            raise ValueError('Nomad response is missing token: \n' + str(response.json()))
        except Exception:
            raise RuntimeError('Something went wrong trying to get Nomad authentication token')

    @with_token
    def create_dataset(self, dataset_name, token=None):
        '''Create a dataset to group a series of NOMAD entries'''
        try:
            response = self.requests.post(
                self.url + 'datasets/',
                headers={'Authorization': f'Bearer {token}', 'Accept': 'application/json'},
                json={"dataset_name": dataset_name},
                timeout=10
            )
            dataset_id = response.json().get('dataset_id')
            if dataset_id:
                return dataset_id

            raise ValueError('Nomad response is missing dataset_id:\n ' + str(response.json()))
        except Exception as e:
            raise RuntimeError('Something went wrong trying to create a Nomad dataset') from e

    @with_token
    def upload(self, upload_file, token=None):
        '''Upload a single file for NOMAD upload, e.g., zip format'''
        def upload():
            try:
                response = self.requests.post(
                    self.url + 'uploads',
                    headers={'Authorization': f'Bearer {token}', 'Accept': 'application/json'},
                    data=f, timeout=30)
                upload_id = response.json().get('upload_id')
                if upload_id:
                    return upload_id
                raise ValueError('Nomad response is missing upload_id:\n ' + str(response.json()))
            except Exception as e:
                raise RuntimeError('Something went wrong uploading to NOMAD') from e
        if isinstance(upload_file, str):
            with open(upload_file, 'rb') as f:
                upload()
        else:
            upload()

    @with_token
    def check_upload_status(self, upload_id, token=None):
        '''
        # upload success => returns 'Process publish_upload completed successfully'
        # publish success => 'Process publish_upload completed successfully'
        '''
        try:
            response = self.requests.get(
                self.url + 'uploads/' + upload_id,
                headers={'Authorization': f'Bearer {token}'}, timeout=30)
            status_message = response.json().get('data').get('last_status_message')
            if status_message:
                return status_message

            raise ValueError('Nomad response is missing status_message: \n' + str(response.json()))
        except Exception as e:
            raise RuntimeError('Something went wrong trying to check the status of upload' + upload_id) from e
            # upload gets deleted from the upload staging area once published...or in this case something went wrong

    @with_token
    def edit_upload_metadata(self, upload_id, metadata, token=None):
        '''
        Example of new metadata:
        ..code-block:: text
            upload_name = 'Test_Upload_Name'
            metadata = {
                "metadata": {
                "upload_name": upload_name,
                "references": ["https://doi.org/xx.xxxx/xxxxxx"],
                "datasets": dataset_id,
                "embargo_length": 0,
                "coauthors": ["coauthor@affiliation.de"],
                "comment": 'This is a test upload...'
                },
            }
        '''

        try:
            response = self.requests.post(
                self.url + 'uploads/' + upload_id + '/edit',
                headers={'Authorization': f'Bearer {token}', 'Accept': 'application/json'},
                json=metadata, timeout=30)
            return response
        except Exception as e:
            raise RuntimeError('Something went wrong trying to add metadata to upload' + upload_id) from e

    def publish_upload(self, upload_id, token=None):
        '''Publish an upload'''
        try:
            response = self.requests.post(
                self.url + 'uploads/' + upload_id + '/action/publish',
                headers={'Authorization': f'Bearer {token}', 'Accept': 'application/json'},
                timeout=30)
            return response
        except Exception as e:
            raise RuntimeError('Something went wrong trying to publish upload: ' + upload_id) from e
