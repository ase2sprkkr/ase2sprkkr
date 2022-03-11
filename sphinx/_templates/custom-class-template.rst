{{ name | escape | underline}}

Full name: **{{ fullname | escape }}**

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
   :members:
   :show-inheritance:
   :inherited-members:

   {% block methods %}
   .. automethod:: __init__

   {% if methods %}
   .. rubric:: {{ _('Methods') }}

   .. autosummary::
      :template: custom-base-template.rst
     {% for item in methods %}
        ~{{ name }}.{{ item }}
     {%- endfor %}

   {% endif %}
   {% endblock %}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: {{ _('Attributes') }}

   .. autosummary::
      :template: custom-base-template.rst
     {% for item in attributes %}
        ~{{ name }}.{{ item }}
     {%- endfor %}

   {% endif %}
   {% endblock %}
