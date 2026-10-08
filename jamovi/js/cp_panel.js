// Panel behaviour shared by every clinpred analysis (jus 3.0 handlers, wired by
// name): the interaction builder's pool, and one row per factor in Reference
// Levels and per linear covariate in Covariate Scaling. Each analysis's
// <basename>.js re-exports this object. jus 3.0 has no context.valuesToItems /
// cloneArray (the jus 2.0 events API): the helpers are the global `utils`.

'use strict';

module.exports = {

    view_loaded: function(ui) {
        this.updateAll(ui);
    },

    view_updated: function(ui) {
        this.updateAll(ui);
    },

    covs_changed: function(ui) {
        this.updateAll(ui);
    },

    splineVars_changed: function(ui) {
        this.updateAll(ui);
    },

    factors_changed: function(ui) {
        this.updateAll(ui);
    },

    refLevels_changed: function(ui) {
        this.updateLevelControls(ui);
    },

    interactionSupplier_updated: function(ui) {
        ui.interactionSupplier.setValue(utils.valuesToItems(this.predictors(ui), FormatDef.variable));
    },

    updateAll: function(ui) {
        this.calcInteractionPool(ui);
        this.updateRefLevels(ui);
        this.updateLevelControls(ui);
        this.updateScales(ui);
        this.updateBaseline(ui);
    },

    predictors: function(ui) {
        let a = utils.clone(ui.covs.value(), []);
        let b = utils.clone(ui.splineVars.value(), []);
        let c = utils.clone(ui.factors.value(), []);
        return a.concat(b).concat(c);
    },

    // terms that reference a removed predictor are pruned
    calcInteractionPool: function(ui) {
        let vars = this.predictors(ui);
        ui.interactionSupplier.setValue(utils.valuesToItems(vars, FormatDef.variable));

        let varsDiff = this.findChanges("predictorList", vars, true, FormatDef.variable);
        let termsList = utils.clone(ui.interactions.value(), []);
        let changed = false;
        for (let i = 0; i < varsDiff.removed.length; i++) {
            for (let j = 0; j < termsList.length; j++) {
                if (FormatDef.term.contains(termsList[j], varsDiff.removed[i])) {
                    termsList.splice(j, 1);
                    changed = true;
                    j -= 1;
                }
            }
        }
        if (changed)
            ui.interactions.setValue(termsList);
    },

    // keep the choice already made for a variable that stays in the list
    syncList: function(list, vars, make) {
        let out = [];
        for (let i = 0; i < vars.length; i++) {
            let found = null;
            for (let j = 0; j < list.length; j++) {
                if (list[j].var === vars[i]) {
                    found = list[j];
                    break;
                }
            }
            out.push(found === null ? make(vars[i]) : found);
        }
        return out;
    },

    updateRefLevels: function(ui) {
        let factors = utils.clone(ui.factors.value(), []);
        let current = utils.clone(ui.refLevels.value(), []);
        ui.refLevels.setValue(this.syncList(current, factors, (v) => ({ var: v, ref: null })));
    },

    // each LevelSelector lists the levels of its own row's variable
    updateLevelControls: function(ui) {
        let list = ui.refLevels.value() || [];
        ui.refLevels.applyToItems(0, (item, index, column) => {
            if (column === 1 && list[index])
                item.setPropertyValue('variable', list[index].var);
        });
    },

    // one row per predictor, with its "baseline" tick kept
    updateBaseline: function(ui) {
        let current = utils.clone(ui.baseline.value(), []);
        ui.baseline.setValue(this.syncList(current, this.predictors(ui), (v) => ({ var: v, forced: false })));
    },

    updateScales: function(ui) {
        let covs = utils.clone(ui.covs.value(), []);
        let current = utils.clone(ui.covScales.value(), []);
        ui.covScales.setValue(this.syncList(current, covs, (v) => ({ var: v, scale: 'unit' })));
    }
};
