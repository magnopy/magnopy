{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. role:: custom-header

.. autoclass:: {{ objname }}

   {% block methods %}
   {% if methods %}
   :custom-header:`{{ _('Methods:') }}`

   .. autosummary::
      :toctree:
   {% for item in methods %}
      {%- if not item.startswith('_') or item in ['__call__'] %}
      ~{{ name }}.{{ item }}
      {%- endif -%}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes %}
   {% if attributes %}
   :custom-header:`{{ _('Properties:') }}`

   .. autosummary::
      :toctree:
   {% for item in attributes %}
      {%- if not item.startswith('_') %}
      ~{{ name }}.{{ item }}
      {%- endif -%}
   {%- endfor %}
   {% endif %}
   {% endblock %}
